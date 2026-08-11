from random import Random

import pysam
import pytest

from generate_demo_data import (
    CNV_LENGTH,
    CNV_SEGMENTS,
    CNV_TUMOUR_PURITY,
    add_sequencing_errors,
    build_cnv_variant_sites,
    write_rna_fusion_bam,
    write_rna_fusion_reference,
)


def test_copy_number_model_links_depth_log2_and_baf_to_state():
    neutral, loss, _, gain, _ = CNV_SEGMENTS

    assert neutral.observed_copy_ratio() == pytest.approx(1.0)
    assert neutral.log2_ratio() == pytest.approx(0.0)
    assert neutral.baf_centres() == pytest.approx((0.5,))

    assert CNV_TUMOUR_PURITY == pytest.approx(0.75)
    assert loss.observed_copy_ratio() == pytest.approx(0.625)
    assert loss.log2_ratio() == pytest.approx(-0.6780719)
    assert loss.baf_centres() == pytest.approx((0.2, 0.8))

    assert gain.observed_copy_ratio() == pytest.approx(1.375)
    assert gain.log2_ratio() == pytest.approx(0.4594316)
    assert gain.baf_centres() == pytest.approx((4 / 11, 7 / 11))


def test_cnv_variant_sites_are_realistically_spaced_and_cover_every_state():
    sites = build_cnv_variant_sites("A" * CNV_LENGTH)

    assert len(sites) == 83
    assert all(left.position < right.position for left, right in zip(sites, sites[1:]))
    assert min(
        sum(site.segment is segment for site in sites)
        for segment in CNV_SEGMENTS
    ) >= 9
    assert all(
        site.target_baf in site.segment.baf_centres()
        for site in sites
    )


def test_error_model_changes_bases_and_marks_them_low_quality():
    sequence = list("ACGTACGT")
    qualities = [36] * len(sequence)

    add_sequencing_errors(sequence, qualities, Random(7), error_rate=1.0)

    assert "".join(sequence) != "ACGTACGT"
    assert all(observed != expected for observed, expected in zip(sequence, "ACGTACGT"))
    assert qualities == [12] * len(sequence)


def test_rna_demo_contains_junction_split_and_spanning_fusion_evidence(tmp_path):
    reference_path = tmp_path / "rna.fa"
    bam_path = tmp_path / "rna.bam"
    references = write_rna_fusion_reference(reference_path)
    write_rna_fusion_bam(bam_path, references)

    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        reads = list(bam.fetch(until_eof=True))

    assert len(reads) == 101
    assert sum(any(op == 3 for op, _ in read.cigartuples or []) for read in reads) == 69
    assert sum(read.has_tag("SA") for read in reads) == 20
    assert sum(read.is_paired for read in reads) == 12
