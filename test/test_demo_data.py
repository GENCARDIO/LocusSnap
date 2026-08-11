from random import Random

import pysam
import pytest

from generate_demo_data import (
    ALK_DNA_BREAKPOINT,
    CNV_LENGTH,
    CNV_SEGMENTS,
    CNV_TUMOUR_PURITY,
    CYP2D6_PHASE_SET,
    CYP2D6_REFERENCE_END,
    CYP2D6_REFERENCE_START,
    CYP2D6_STAR4_VARIANTS,
    EML4_DNA_BREAKPOINT,
    add_sequencing_errors,
    build_cnv_variant_sites,
    write_eml4_alk_dna_bam,
    write_cyp2d6_bam,
    write_cyp2d6_vcf,
    write_rna_fusion_bam,
    write_rna_fusion_reference,
)
from locus_snap.read_model import compute_pair_orientation


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


def test_eml4_alk_demo_contains_rna_fusion_and_dna_inversion_evidence(tmp_path):
    reference_path = tmp_path / "rna.fa"
    rna_bam_path = tmp_path / "rna.bam"
    dna_bam_path = tmp_path / "dna.bam"
    references = write_rna_fusion_reference(reference_path)
    write_rna_fusion_bam(rna_bam_path, references)
    write_eml4_alk_dna_bam(dna_bam_path, references)

    with pysam.AlignmentFile(str(rna_bam_path), "rb") as bam:
        rna_reads = list(bam.fetch(until_eof=True))

    assert len(rna_reads) == 150
    assert sum(any(op == 3 for op, _ in read.cigartuples or []) for read in rna_reads) == 94
    assert sum(read.has_tag("SA") for read in rna_reads) == 36
    assert sum(read.is_paired for read in rna_reads) == 20
    assert all(read.reference_name == "chr2" for read in rna_reads)
    assert sum("EML4e13_ALKe20" in read.query_name for read in rna_reads) == 56

    with pysam.AlignmentFile(str(dna_bam_path), "rb") as bam:
        dna_reads = list(bam.fetch(until_eof=True))

    assert len(dna_reads) == 582
    assert sum(read.has_tag("SA") for read in dna_reads) == 20
    orientations = {
        compute_pair_orientation(read)
        for read in dna_reads
        if read.is_paired and not read.is_proper_pair
    }
    assert {"FF", "RR"} <= orientations
    assert EML4_DNA_BREAKPOINT - ALK_DNA_BREAKPOINT == 13_075_342


def test_cyp2d6_demo_links_phased_reads_to_star4_vcf(tmp_path):
    bam_path = tmp_path / "cyp2d6.bam"
    vcf_path = tmp_path / "cyp2d6.vcf"
    reference = ["A"] * (CYP2D6_REFERENCE_END - CYP2D6_REFERENCE_START)
    for position, ref, *_ in CYP2D6_STAR4_VARIANTS:
        reference[position - 1 - CYP2D6_REFERENCE_START] = ref

    write_cyp2d6_bam(bam_path, "".join(reference))
    write_cyp2d6_vcf(vcf_path, bam_path)

    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        reads = list(bam.fetch("chr22"))
    assert reads
    assert {read.get_tag("HP") for read in reads} == {1, 2}
    assert {read.get_tag("PS") for read in reads} == {CYP2D6_PHASE_SET}

    with pysam.VariantFile(str(vcf_path)) as vcf:
        records = list(vcf)
    assert len(records) == len(CYP2D6_STAR4_VARIANTS) == 5
    assert all(record.samples["PGX_DEMO"]["GT"] == (1, 0) for record in records)
    assert all(record.samples["PGX_DEMO"].phased for record in records)
    assert all(record.samples["PGX_DEMO"]["DP"] > 0 for record in records)
    assert all(min(record.samples["PGX_DEMO"]["AD"]) > 0 for record in records)
    assert sum(record.info["ROLE"] == "core" for record in records) == 1
