import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from locus_snap.snapshot import render_multi_locus_snapshots
from locus_snap.cli import main


TEST_BAM = os.path.join(os.path.dirname(__file__), "test.bam")


def test_explicit_loci_render_independent_scales_and_stacked_samples(tmp_path):
    output, summary = render_multi_locus_snapshots(
        bam_paths=[TEST_BAM, TEST_BAM],
        regions=[
            ("chr9", 101867480, 101867550),
            ("chr9", 101867550, 101867620),
        ],
        sample_labels=["Tumour", "Normal"],
        region_labels=["Left breakpoint", "Right breakpoint"],
        output_dir=str(tmp_path),
        output_name="explicit-loci.svg",
        show_ideogram=False,
        show_coverage=False,
        show_legend=False,
        link_breakpoints=True,
        max_rows=1,
        fig_width=7,
        dpi=40,
    )

    svg = (tmp_path / "explicit-loci.svg").read_text(encoding="utf-8")
    assert output.endswith("explicit-loci.svg")
    assert "Left breakpoint" in svg
    assert "Right breakpoint" in svg
    assert "chr9:101,867,481-101,867,550" in svg
    assert "chr9:101,867,551-101,867,620" in svg
    assert "Tumour" in svg
    assert "Normal" in svg
    assert "#7a1f5c" in svg
    assert "Tumour · Left breakpoint" in summary
    assert "Normal · Right breakpoint" in summary


def test_explicit_loci_support_three_panels(tmp_path):
    output, _summary = render_multi_locus_snapshots(
        bam_paths=[TEST_BAM],
        regions=[
            ("chr9", 101867480, 101867520),
            ("chr9", 101867520, 101867560),
            ("chr9", 101867560, 101867600),
        ],
        output_dir=str(tmp_path),
        output_name="three-loci.png",
        show_ideogram=False,
        show_coverage=False,
        show_legend=False,
        link_breakpoints=True,
        max_rows=1,
        fig_width=8,
        dpi=40,
    )

    assert output.endswith("three-loci.png")
    assert os.path.getsize(output) > 0


def test_breakpoint_links_require_repeated_regions(caplog):
    return_code = main([
        "--bam", TEST_BAM,
        "--region", "chr9:101867481-101867520",
        "--link_breakpoints",
        "--refseq", "none",
    ])

    assert return_code == 1
    assert "requires at least two --region" in caplog.text
