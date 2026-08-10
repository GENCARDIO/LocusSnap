import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from locus_snap.cli import main, parse_tag_colors


TEST_BAM = os.path.join(os.path.dirname(__file__), "test.bam")


def test_parse_tag_color_overrides_validates_and_normalizes_values():
    assert parse_tag_colors([
        "tumour=#445566", "UNTAGGED=gold",
    ]) == {"tumour": "#445566", "untagged": "gold"}
    with pytest.raises(ValueError, match="VALUE=COLOR"):
        parse_tag_colors(["tumour"])
    with pytest.raises(ValueError, match="valid Matplotlib colour"):
        parse_tag_colors(["tumour=definitely-not-a-colour"])


def test_cli_groups_and_colours_reads_by_generic_tag(tmp_path):
    rc = main([
        "--bam", TEST_BAM,
        "--region", "chr9:101867481-101867620",
        "--group_by_tag", "RG",
        "--tag_filter", "untagged",
        "--tag_label", "Library",
        "--tag_color", "untagged=#777777",
        "--display_mode", "collapse",
        "--no_coverage",
        "--no_ideogram",
        "--refseq", "none",
        "--output_dir", str(tmp_path),
        "--output_name", "tag-groups.svg",
        "--fig_width", "7",
        "--dpi", "40",
    ])

    assert rc == 0
    svg = (tmp_path / "tag-groups.svg").read_text(encoding="utf-8")
    assert "Library=untagged" in svg
    assert "#777777" in svg


def test_cli_rejects_tag_controls_without_tag_mode(tmp_path):
    rc = main([
        "--bam", TEST_BAM,
        "--region", "chr9:101867481-101867620",
        "--tag_filter", "untagged",
        "--refseq", "none",
        "--output_dir", str(tmp_path),
    ])
    assert rc == 1


def test_cli_rejects_generic_tag_and_haplotype_colour_modes(tmp_path):
    rc = main([
        "--bam", TEST_BAM,
        "--region", "chr9:101867481-101867620",
        "--color_by_tag", "RG",
        "--haplotype_view", "color",
        "--refseq", "none",
        "--output_dir", str(tmp_path),
    ])
    assert rc == 1
