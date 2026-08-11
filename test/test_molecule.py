from array import array
import os
import sys

import pysam
import pytest
from matplotlib.colors import to_hex

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from locus_snap.config import DEFAULT_MOLECULE_COLORS
from locus_snap.molecule import (
    build_molecule_consensus_reads,
    resolve_molecule_tag,
)
from locus_snap.read_model import AlignedRead
from locus_snap.render import AlignmentRenderer, compute_coverage, compute_snv_evidence
from locus_snap.snapshot import BamSnapshot


HEADER = pysam.AlignmentHeader.from_references(["chr1"], [10_000])


class FakeReference:
    def __init__(self, start=0, sequence="A" * 10_000):
        self.start = start
        self.sequence = sequence
        self.available = True

    def base_at(self, position):
        index = position - self.start
        if 0 <= index < len(self.sequence):
            return self.sequence[index]
        return None


def make_read(
    name, sequence="AAAAAAAAAA", start=100, *, umi="UMI-1", cell=None,
    molecule_id=None, duplicate=False, reverse=False, read_side=0,
    quality=35, reference=None,
):
    segment = pysam.AlignedSegment(HEADER)
    segment.query_name = name
    segment.query_sequence = sequence
    flag = 0
    if duplicate:
        flag |= 0x400
    if reverse:
        flag |= 0x10
    if read_side:
        flag |= 0x1 | (0x40 if read_side == 1 else 0x80)
    segment.flag = flag
    segment.reference_id = 0
    segment.reference_start = start
    segment.mapping_quality = 60
    segment.cigartuples = [(0, len(sequence))]
    segment.query_qualities = array("B", [quality] * len(sequence))
    if umi is not None:
        segment.set_tag("RX", umi)
    if cell is not None:
        segment.set_tag("CB", cell)
    if molecule_id is not None:
        segment.set_tag("MI", molecule_id)
    return AlignedRead(segment, reference=reference)


def test_auto_tag_prefers_mi_on_equal_coverage_and_validates_explicit_tag():
    reads = [
        make_read("one", molecule_id="mol-1"),
        make_read("two", molecule_id="mol-2"),
    ]

    assert resolve_molecule_tag(reads) == "MI"
    assert resolve_molecule_tag(reads, "rx") == "RX"
    with pytest.raises(ValueError, match="molecule tag UB"):
        resolve_molecule_tag(reads, "UB")
    with pytest.raises(ValueError, match="auto or one of"):
        resolve_molecule_tag(reads, "XX")


def test_positional_umi_families_respect_cell_barcode_read_side_and_tolerance():
    reads = [
        make_read("a", start=100, cell="CELL-A", read_side=1),
        make_read("b", start=101, cell="CELL-A", read_side=1),
        make_read("other-cell", start=100, cell="CELL-B", read_side=1),
        make_read("other-side", start=100, cell="CELL-A", read_side=2),
        make_read("far", start=106, cell="CELL-A", read_side=1),
    ]

    result = build_molecule_consensus_reads(
        reads, requested_tag="RX", position_tolerance=2
    )

    assert result.family_count == 4
    assert sorted(read.molecule_family_size for read in result.reads) == [1, 1, 1, 2]
    assert any("CELL-A:UMI-1" in read.molecule_id for read in result.reads)


def test_consensus_uses_majority_base_and_reports_family_metadata():
    reference = FakeReference()
    sequences = ["AAAAGAAAAA", "AAAAGAAAAA", "AAAAAAAAAA"]
    reads = [
        make_read(
            f"read-{index}", sequence=sequence, duplicate=index > 0,
            reverse=index == 2, reference=reference,
        )
        for index, sequence in enumerate(sequences)
    ]

    result = build_molecule_consensus_reads(
        reads, requested_tag="RX", minimum_family_size=2,
        minimum_consensus_fraction=0.60, reference=reference,
    )
    consensus = result.reads[0]

    assert result.input_read_count == 3
    assert result.tagged_read_count == 3
    assert result.family_count == 1
    assert result.filtered_family_count == 0
    assert consensus.molecule_family_size == 3
    assert consensus.molecule_duplicate_reads == 2
    assert consensus.molecule_is_duplex is True
    assert consensus.molecule_member_names == ("read-0", "read-1", "read-2")
    assert consensus.base_at(104) == "G"
    assert consensus.mismatches == [(104, "G")]
    assert consensus.molecule_consensus_fraction == pytest.approx(29 / 30)


def test_minimum_family_size_filters_singletons():
    reads = [
        make_read("family-a1", umi="A"),
        make_read("family-a2", umi="A"),
        make_read("singleton", umi="B", start=200),
    ]

    result = build_molecule_consensus_reads(
        reads, requested_tag="RX", minimum_family_size=2
    )

    assert len(result.reads) == 1
    assert result.reads[0].molecule_family_size == 2
    assert result.family_count == 2
    assert result.filtered_family_count == 1


def test_snapshot_retains_duplicate_members_inside_molecule_families(tmp_path):
    bam_path = tmp_path / "molecules.bam"
    reads = [
        make_read("original", umi="FAMILY"),
        make_read("duplicate-1", umi="FAMILY", duplicate=True),
        make_read("duplicate-2", umi="FAMILY", duplicate=True),
    ]
    with pysam.AlignmentFile(str(bam_path), "wb", header=HEADER) as bam:
        for read in reads:
            bam.write(read.segment)
    pysam.index(str(bam_path))

    snapshot = BamSnapshot(
        str(bam_path), "chr1", 95, 120,
        output_dir=str(tmp_path), molecule_mode=True,
        molecule_tag="RX", min_family_size=2,
    )
    consensus_reads = snapshot.load_reads()

    assert len(consensus_reads) == 1
    assert consensus_reads[0].molecule_family_size == 3
    assert consensus_reads[0].molecule_duplicate_reads == 2


def test_snapshot_rejects_incompatible_molecule_views(tmp_path):
    common = dict(
        bam="unused.bam", chrom="chr1", start=100, end=110,
        output_dir=str(tmp_path), molecule_mode=True,
    )
    with pytest.raises(ValueError, match="Molecule mode"):
        BamSnapshot(**common, view_as_pairs=True)
    with pytest.raises(ValueError, match="Molecule mode"):
        BamSnapshot(**common, show_base_modifications=True)
    with pytest.raises(ValueError, match="Molecule mode"):
        BamSnapshot(**common, haplotype_view="color")
    with pytest.raises(ValueError, match="Molecule mode"):
        BamSnapshot(**common, read_tag="RG", tag_view="split")


def test_consensus_family_counts_once_for_coverage_and_vaf():
    reference = FakeReference()
    reads = [
        make_read("a", sequence="AAAAGAAAAA", reference=reference),
        make_read("b", sequence="AAAAGAAAAA", duplicate=True, reference=reference),
        make_read("c", sequence="AAAAGAAAAA", duplicate=True, reference=reference),
    ]
    consensus_reads = build_molecule_consensus_reads(
        reads, requested_tag="RX", reference=reference
    ).reads

    depths = compute_coverage(consensus_reads, 100, 110)
    snv_depth, evidence = compute_snv_evidence(consensus_reads, 100, 110)

    assert depths == [1] * 10
    assert snv_depth == [1] * 10
    assert evidence[104]["G"].count == 1


def test_renderer_colors_singleton_consensus_and_duplex_molecules():
    renderer = AlignmentRenderer(molecule_mode=True, shade_by_mapq=False)
    singleton = make_read("singleton", umi="A")
    family = make_read("family", umi="B")
    family.molecule_family_size = 3
    duplex = make_read("duplex", umi="C")
    duplex.molecule_family_size = 4
    duplex.molecule_is_duplex = True

    assert to_hex(renderer.read_style(singleton)[0]) == DEFAULT_MOLECULE_COLORS["singleton"]
    assert to_hex(renderer.read_style(family)[0]) == DEFAULT_MOLECULE_COLORS["consensus"]
    assert to_hex(renderer.read_style(duplex)[0]) == DEFAULT_MOLECULE_COLORS["duplex"]


@pytest.mark.parametrize(
    "kwargs, message",
    [
        ({"minimum_family_size": 0}, "at least one"),
        ({"position_tolerance": -1}, "cannot be negative"),
        ({"minimum_consensus_fraction": 0.49}, "between 0.5 and 1"),
    ],
)
def test_molecule_configuration_validation(kwargs, message):
    with pytest.raises(ValueError, match=message):
        build_molecule_consensus_reads([make_read("one")], **kwargs)
