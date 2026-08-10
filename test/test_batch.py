import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from locus_snap.batch import (
    BatchRegion, BatchResult, looks_like_vcf, parse_bed_regions, parse_vcf_regions,
    write_html_report,
)
from locus_snap.cli import main
from locus_snap.metrics import RegionSummary

TEST_BAM = os.path.join(os.path.dirname(__file__), "test.bam")
VCF_HEADER = (
    "##fileformat=VCFv4.2\n"
    '##INFO=<ID=END,Number=1,Type=Integer,Description="End position">\n'
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
)


def write_bed(tmp_path, text, name="regions.bed"):
    path = tmp_path / name
    path.write_text(text, encoding="utf-8")
    return str(path)


def write_vcf(tmp_path, records_text, name="regions.vcf"):
    path = tmp_path / name
    path.write_text(VCF_HEADER + records_text, encoding="utf-8")
    return str(path)


def test_parse_bed_regions_bed3_and_bed4(tmp_path):
    bed = write_bed(
        tmp_path,
        "# comment\n"
        "track name=demo\n"
        "\n"
        "chr1\t100\t200\n"
        "chr2\t300\t400\tmy_region\n",
    )
    regions = parse_bed_regions(bed)

    assert regions == [
        BatchRegion(chrom="chr1", start=100, end=200, name="chr1_100_200"),
        BatchRegion(chrom="chr2", start=300, end=400, name="my_region"),
    ]
    assert regions[1].display == "chr2:301-400"


def test_parse_bed_regions_applies_flank(tmp_path):
    bed = write_bed(tmp_path, "chr1\t100\t200\n")
    regions = parse_bed_regions(bed, flank=50)
    assert regions == [BatchRegion(chrom="chr1", start=50, end=250, name="chr1_50_250")]


def test_parse_bed_regions_flank_does_not_go_negative(tmp_path):
    bed = write_bed(tmp_path, "chr1\t10\t20\n")
    regions = parse_bed_regions(bed, flank=1000)
    assert regions[0].start == 0


def test_parse_bed_regions_deduplicates_names(tmp_path):
    bed = write_bed(tmp_path, "chr1\t100\t200\tdup\nchr1\t300\t400\tdup\n")
    regions = parse_bed_regions(bed)
    assert [r.name for r in regions] == ["dup", "dup_2"]


def test_parse_bed_regions_rejects_bad_lines(tmp_path):
    bed = write_bed(tmp_path, "chr1\t100\n")
    with pytest.raises(ValueError, match="expected at least 3 BED columns"):
        parse_bed_regions(bed)

    bed = write_bed(tmp_path, "chr1\tabc\t200\n", name="bad_int.bed")
    with pytest.raises(ValueError, match="must be integers"):
        parse_bed_regions(bed)

    bed = write_bed(tmp_path, "chr1\t200\t100\n", name="bad_order.bed")
    with pytest.raises(ValueError, match="must be greater than start"):
        parse_bed_regions(bed)


def test_parse_bed_regions_rejects_empty_file(tmp_path):
    bed = write_bed(tmp_path, "# only comments\n\n")
    with pytest.raises(ValueError, match="no regions found"):
        parse_bed_regions(bed)


def test_parse_bed_regions_missing_file():
    with pytest.raises(ValueError, match="Cannot find --batch_regions file"):
        parse_bed_regions("/no/such/file.bed")


def _summary(n_reads=10, n_gapped=2, n_discordant=1, n_softclipped=3):
    return RegionSummary(
        label="r", n_reads=n_reads, n_gapped=n_gapped, n_long_gap=0, long_gap_threshold=10,
        max_gap=5, mean_gap_of_gapped=5.0, total_gap_bp=10, n_with_sa=0, n_cross_chrom_sa=0,
        n_discordant=n_discordant, n_interchrom=0, n_softclipped=n_softclipped, mean_mapq=50.0,
    )


def test_write_html_report_embeds_images_and_lists_failures(tmp_path):
    image_path = tmp_path / "region_a.png"
    image_path.write_bytes(b"\x89PNG\r\n\x1a\nfake-png-bytes")

    results = [
        BatchResult(
            region=BatchRegion(chrom="chr1", start=0, end=100, name="region_a"),
            output_path=str(image_path), summary=_summary(),
        ),
        BatchResult(
            region=BatchRegion(chrom="chr2", start=0, end=100, name="region_b"),
            error="invalid contig `chr2`",
        ),
    ]
    report_path = tmp_path / "report.html"

    write_html_report(results, str(report_path), "png")
    content = report_path.read_text(encoding="utf-8")

    assert "data:image/png;base64," in content
    assert "region_a" in content
    assert "region_b" in content
    assert "invalid contig" in content
    assert "1 rendered, 1 failed" in content


def test_write_html_report_groups_multiple_samples_and_shows_deltas(tmp_path):
    image_path = tmp_path / "comparison.png"
    image_path.write_bytes(b"\x89PNG\r\n\x1a\nfake-png-bytes")
    tumour = _summary(n_reads=12, n_gapped=3, n_discordant=4, n_softclipped=5)
    normal = _summary(n_reads=8, n_gapped=1, n_discordant=1, n_softclipped=2)
    result = BatchResult(
        region=BatchRegion(chrom="chr1", start=99, end=100, name="variant_a"),
        output_path=str(image_path), summaries=[tumour, normal], summary=tumour,
    )
    report_path = tmp_path / "multi.html"

    write_html_report(
        [result], str(report_path), "png", sample_labels=["Tumour", "Normal"]
    )
    content = report_path.read_text(encoding="utf-8")

    assert "colspan='4'>Tumour" in content
    assert "colspan='4'>Normal" in content
    assert "2 samples" in content
    assert "&Delta; reads" in content
    assert "-4" in content
    assert "baseline" in content


def test_write_html_report_rejects_non_browser_format(tmp_path):
    results = [BatchResult(region=BatchRegion(chrom="chr1", start=0, end=100, name="r"))]
    with pytest.raises(ValueError, match="browser-viewable"):
        write_html_report(results, str(tmp_path / "report.html"), "pdf")


def test_batch_cli_renders_all_regions_and_report(tmp_path):
    bed = write_bed(
        tmp_path,
        "chr9\t101867480\t101867620\tsite_a\n"
        "chr9\t101867500\t101867560\n"
        "chrZZZ\t1\t100\tbad_contig\n",
    )
    output_dir = tmp_path / "out"

    rc = main([
        "--bam", TEST_BAM,
        "--batch_regions", bed,
        "--report",
        "--refseq", "none",
        "--output_dir", str(output_dir),
    ])

    assert rc == 1  # one region failed
    assert (output_dir / "site_a.png").is_file()
    assert (output_dir / "chr9_101867500_101867560.png").is_file()
    assert not (output_dir / "bad_contig.png").exists()
    report = (output_dir / "report.html").read_text(encoding="utf-8")
    assert "site_a" in report
    assert "bad_contig" in report
    assert "data:image/png;base64," in report


def test_batch_cli_all_regions_succeed_returns_zero(tmp_path):
    bed = write_bed(tmp_path, "chr9\t101867480\t101867620\tsite_a\n")
    output_dir = tmp_path / "out"

    rc = main([
        "--bam", TEST_BAM,
        "--batch_regions", bed,
        "--refseq", "none",
        "--output_dir", str(output_dir),
    ])

    assert rc == 0
    assert (output_dir / "site_a.png").is_file()
    assert not (output_dir / "report.html").exists()


def test_report_without_batch_regions_is_rejected(tmp_path):
    rc = main([
        "--bam", TEST_BAM,
        "--region", "chr9:101867481-101867620",
        "--report",
        "--refseq", "none",
        "--output_dir", str(tmp_path),
    ])
    assert rc == 1


def test_batch_regions_renders_multiple_bams_in_parallel_and_reports_samples(tmp_path):
    bed = write_bed(
        tmp_path,
        "chr9\t101867480\t101867620\tsite_a\n"
        "chrZZZ\t1\t100\tbad_contig\n"
        "chr9\t101867500\t101867560\tsite_b\n",
    )
    output_dir = tmp_path / "out"
    rc = main([
        "--bam", TEST_BAM,
        "--bam", TEST_BAM,
        "--sample_label", "Tumour",
        "--sample_label", "Normal",
        "--batch_regions", bed,
        "--report",
        "--threads", "2",
        "--max_rows", "2",
        "--no_coverage",
        "--no_legend",
        "--no_ideogram",
        "--refseq", "none",
        "--output_dir", str(output_dir),
    ])

    assert rc == 1
    assert (output_dir / "site_a.png").is_file()
    assert (output_dir / "site_b.png").is_file()
    assert not (output_dir / "bad_contig.png").exists()
    report = (output_dir / "report.html").read_text(encoding="utf-8")
    assert "colspan='4'>Tumour" in report
    assert "colspan='4'>Normal" in report
    assert report.index("site_a") < report.index("bad_contig") < report.index("site_b")


def test_batch_regions_rejects_non_positive_threads(tmp_path):
    bed = write_bed(tmp_path, "chr9\t101867480\t101867620\tsite_a\n")
    rc = main([
        "--bam", TEST_BAM,
        "--batch_regions", bed,
        "--threads", "0",
        "--refseq", "none",
        "--output_dir", str(tmp_path),
    ])
    assert rc == 1


def test_batch_regions_rejects_sort_base_position(tmp_path):
    bed = write_bed(tmp_path, "chr9\t101867480\t101867620\tsite_a\n")
    rc = main([
        "--bam", TEST_BAM,
        "--batch_regions", bed,
        "--sort_by", "base",
        "--sort_base_position", "101867500",
        "--refseq", "none",
        "--output_dir", str(tmp_path),
    ])
    assert rc == 1


def test_batch_regions_rejects_metrics_tsv(tmp_path):
    bed = write_bed(tmp_path, "chr9\t101867480\t101867620\tsite_a\n")
    rc = main([
        "--bam", TEST_BAM,
        "--batch_regions", bed,
        "--metrics_tsv", str(tmp_path / "out.tsv"),
        "--refseq", "none",
        "--output_dir", str(tmp_path),
    ])
    assert rc == 1


def test_looks_like_vcf_dispatch():
    assert looks_like_vcf("calls.vcf")
    assert looks_like_vcf("calls.vcf.gz")
    assert looks_like_vcf("calls.VCF.GZ")
    assert looks_like_vcf("calls.bcf")
    assert not looks_like_vcf("regions.bed")
    assert not looks_like_vcf("regions.txt")


def test_parse_vcf_regions_uses_id_column_when_present(tmp_path):
    vcf = write_vcf(tmp_path, "chr1\t100\trs123\tA\tG\t.\tPASS\t.\n")
    regions = parse_vcf_regions(vcf)
    assert regions == [BatchRegion(chrom="chr1", start=99, end=100, name="rs123")]


def test_parse_vcf_regions_falls_back_to_chrom_pos_ref_alt_name(tmp_path):
    vcf = write_vcf(tmp_path, "chr1\t100\t.\tA\tG\t.\tPASS\t.\n")
    regions = parse_vcf_regions(vcf)
    assert regions == [BatchRegion(chrom="chr1", start=99, end=100, name="chr1_100_A_G")]


def test_parse_vcf_regions_multiallelic_line_is_one_region(tmp_path):
    vcf = write_vcf(tmp_path, "chr1\t100\t.\tAT\tA,ATT\t.\tPASS\t.\n")
    regions = parse_vcf_regions(vcf)
    assert len(regions) == 1
    assert regions[0].name == "chr1_100_AT_A_ATT"  # comma sanitized to '_'


def test_parse_vcf_regions_resolves_symbolic_end_from_info(tmp_path):
    vcf = write_vcf(tmp_path, "chr1\t100\t.\tN\t<DEL>\t.\tPASS\tEND=200\n")
    regions = parse_vcf_regions(vcf)
    # '<' and '>' are sanitized to '_' independently
    assert regions == [BatchRegion(chrom="chr1", start=99, end=200, name="chr1_100_N__DEL_")]


def test_parse_vcf_regions_skips_records_without_alt(tmp_path):
    vcf = write_vcf(
        tmp_path,
        "chr1\t100\trs1\tA\tG\t.\tPASS\t.\n"
        "chr1\t200\trs2\tA\t.\t.\tPASS\t.\n",
    )
    regions = parse_vcf_regions(vcf)
    assert [r.name for r in regions] == ["rs1"]


def test_parse_vcf_regions_applies_flank(tmp_path):
    vcf = write_vcf(tmp_path, "chr1\t100\trs1\tA\tG\t.\tPASS\t.\n")
    regions = parse_vcf_regions(vcf, flank=10)
    assert regions == [BatchRegion(chrom="chr1", start=89, end=110, name="rs1")]


def test_parse_vcf_regions_rejects_empty_file(tmp_path):
    vcf = write_vcf(tmp_path, "")
    with pytest.raises(ValueError, match="no variants with an ALT allele"):
        parse_vcf_regions(vcf)


def test_parse_vcf_regions_missing_file():
    with pytest.raises(ValueError, match="Cannot find --batch_regions file"):
        parse_vcf_regions("/no/such/file.vcf")


def test_batch_cli_renders_from_vcf_and_report(tmp_path):
    vcf = write_vcf(
        tmp_path,
        "chr9\t101867500\trs_snv\tA\tG\t.\tPASS\t.\n"
        "chr9\t101867560\t.\tN\t<DEL>\t.\tPASS\tEND=101867600\n",
    )
    output_dir = tmp_path / "out"

    rc = main([
        "--bam", TEST_BAM,
        "--batch_regions", vcf,
        "--flank", "20",
        "--report",
        "--refseq", "none",
        "--output_dir", str(output_dir),
    ])

    assert rc == 0
    assert (output_dir / "rs_snv.png").is_file()
    assert (output_dir / "chr9_101867560_N__DEL_.png").is_file()
    report = (output_dir / "report.html").read_text(encoding="utf-8")
    assert "rs_snv" in report
    assert "data:image/png;base64," in report


def test_batch_cli_still_dispatches_bed_files_to_bed_parser(tmp_path):
    bed = write_bed(tmp_path, "chr9\t101867480\t101867620\tsite_a\n")
    output_dir = tmp_path / "out"

    rc = main([
        "--bam", TEST_BAM,
        "--batch_regions", bed,
        "--refseq", "none",
        "--output_dir", str(output_dir),
    ])

    assert rc == 0
    assert (output_dir / "site_a.png").is_file()
