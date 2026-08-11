import os
import sys
from types import SimpleNamespace

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from locus_snap.annotations import AnnotationItem, LoadedAnnotationTrack
from locus_snap.read_model import CigarBlock, SAEntry
from locus_snap.rna import (
    annotated_junctions,
    collect_fusion_evidence,
    collect_splice_junctions,
    gene_label_at,
    is_annotated_junction,
    rna_evidence_rows,
    splice_motif,
    transcript_strand,
    write_rna_evidence_tsv,
)


def junction_read(name, strand="+", left_anchor=20, right_anchor=25, read2=False):
    return SimpleNamespace(
        query_name=name, strand=strand, is_read2=read2,
        blocks=[
            CigarBlock("M", 100, 0, left_anchor),
            CigarBlock("N", 120, left_anchor, 80),
            CigarBlock("M", 200, left_anchor, right_anchor),
        ],
        deletions=[(120, 80, True)],
    )


def fusion_read(
    name, *, start=100, end=150, sa_start=500, sa_end=550,
    sa_chrom="chr2", duplicate_record=False, paired=False,
):
    return SimpleNamespace(
        query_name=name, ref_start=start, ref_end=end, strand="+",
        is_reverse=False, soft_clip_left=0, soft_clip_right=50,
        hard_clip_left=0, hard_clip_right=0, mapq=60,
        sa_entries=[SAEntry(sa_chrom, sa_start, sa_end, "+", "50S50M", 55, 1)],
        mate_chrom="chr3" if paired else None,
        mate_start=900 if paired else None,
        mate_is_unmapped=not paired, is_read2=duplicate_record,
    )


class FakeReference:
    available = True

    def __init__(self, bases):
        self.bases = bases

    def base_at(self, position):
        return self.bases.get(position, "A")


def test_junction_collection_deduplicates_names_and_enforces_both_anchors():
    reads = [
        junction_read("same-read"),
        junction_read("same-read"),
        junction_read("short", left_anchor=4),
        junction_read("other", strand="-"),
    ]

    combined = collect_splice_junctions(reads, minimum_anchor=8)
    split = collect_splice_junctions(
        reads, strand_mode="split", minimum_anchor=8
    )

    assert combined[(120, 200, ".")].count == 2
    assert split[(120, 200, "+")].read_names == ("same-read",)
    assert split[(120, 200, "-")].read_names == ("other",)


def test_library_strandness_normalizes_paired_read_orientation():
    read1 = junction_read("r1", strand="+", read2=False)
    read2 = junction_read("r2", strand="-", read2=True)

    assert transcript_strand(read1, "forward") == "+"
    assert transcript_strand(read2, "forward") == "+"
    assert transcript_strand(read1, "reverse") == "-"
    assert transcript_strand(read2, "reverse") == "-"
    assert transcript_strand(read1, "unstranded") == "."


def test_visible_gene_models_classify_known_and_novel_junctions():
    item = AnnotationItem(
        90, 230, "TX1", "+", blocks=[(90, 120), (200, 230)],
        group_label="GENE1",
    )
    track = LoadedAnnotationTrack(
        "Genes", "gtf", "#123456", [item], [[item]], chrom="chr1"
    )
    known = annotated_junctions([track])
    evidence = collect_splice_junctions([junction_read("known")])[(120, 200, ".")]

    assert ("chr1", 120, 200, "+") in known
    assert is_annotated_junction("chr1", evidence, known)
    assert gene_label_at([track], "chr1", 119) == "GENE1"
    assert gene_label_at([track], "chr2", 119) is None


def test_splice_motif_is_reported_in_transcript_orientation():
    plus = collect_splice_junctions(
        [junction_read("plus")], strand_mode="split"
    )[(120, 200, "+")]
    minus = collect_splice_junctions(
        [junction_read("minus", strand="-")], strand_mode="split"
    )[(120, 200, "-")]

    assert splice_motif(
        plus, FakeReference({120: "G", 121: "T", 198: "A", 199: "G"})
    ) == "GT-AG"
    assert splice_motif(
        minus, FakeReference({120: "C", 121: "T", 198: "A", 199: "C"})
    ) == "GT-AG"


def test_fusion_evidence_clusters_jitter_and_deduplicates_split_records():
    reads = [
        fusion_read("split-1", sa_start=500, sa_end=550),
        fusion_read("split-2", start=101, end=151, sa_start=503, sa_end=553),
        fusion_read("split-1", sa_start=500, sa_end=550, duplicate_record=True),
        fusion_read("split-and-pair", paired=True),
    ]

    evidence = collect_fusion_evidence(
        reads, "chr1", breakpoint_tolerance=5, minimum_distance=100
    )

    chr2 = next(item for item in evidence if item.partner_chrom == "chr2")
    chr3 = next(item for item in evidence if item.partner_chrom == "chr3")
    assert chr2.local_breakpoint == 150
    assert chr2.partner_breakpoint == 500
    assert chr2.split_reads == 3
    assert chr2.support == 3
    assert chr3.spanning_pairs == 1
    assert chr3.support == 1


def test_rna_evidence_rows_and_tsv_export_junctions_and_fusions(tmp_path):
    transcript = AnnotationItem(
        90, 230, "TX1", "+", blocks=[(90, 120), (200, 230)]
    )
    track = LoadedAnnotationTrack(
        "Genes", "gtf", "#123456", [transcript], [[transcript]], chrom="chr1"
    )
    reads = [junction_read("junction"), fusion_read("fusion")]
    rows = rna_evidence_rows(
        reads, "chr1", genomic_tracks=[track], minimum_fusion_reads=1,
        fusion_minimum_distance=100,
    )

    junction = next(row for row in rows if row["event_type"] == "junction")
    fusion = next(row for row in rows if row["event_type"] == "fusion")
    assert (junction["status"], junction["support"]) == ("annotated", 1)
    assert (fusion["partner_chrom"], fusion["split_reads"]) == ("chr2", 1)

    path = tmp_path / "rna.tsv"
    write_rna_evidence_tsv(
        str(path), reads, "chr1", genomic_tracks=[track],
        minimum_fusion_reads=1, fusion_minimum_distance=100,
    )
    text = path.read_text(encoding="utf-8")
    assert text.startswith("event_type\tchrom\tstart\tend")
    assert "junction\tchr1" in text
    assert "fusion\tchr1" in text


@pytest.mark.parametrize(
    "kwargs, message",
    [
        ({"breakpoint_tolerance": -1}, "tolerance"),
        ({"minimum_distance": -1}, "distance"),
        ({"minimum_mapq": -1}, "MAPQ"),
    ],
)
def test_fusion_configuration_validation(kwargs, message):
    with pytest.raises(ValueError, match=message):
        collect_fusion_evidence([], "chr1", **kwargs)
