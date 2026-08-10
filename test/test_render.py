import os
import sys
from types import SimpleNamespace

import matplotlib.pyplot as plt
import pytest
from matplotlib.colors import to_hex
from matplotlib.patches import Polygon

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from locus_snap.annotations import AnnotationItem, LoadedAnnotationTrack
from locus_snap.config import DEFAULT_BASE_COLORS, DEFAULT_VISUAL_COLORS, load_config
from locus_snap.cytobands import Cytoband
from locus_snap.read_model import CigarBlock
from locus_snap.render import (
    AlignmentRenderer,
    HighlightRegion,
    compute_coverage,
    compute_binned_coverage,
    compute_feature_density,
    compute_sparse_snv_evidence,
    compute_snv_evidence,
    compute_snv_counts,
    compute_splice_junctions,
    ellipsize,
    format_scale_length,
    genomic_tick_labels,
    haplotype_color,
    nice_scale_length,
    tag_color,
)


def test_feature_density_counts_interval_overlap_per_bin():
    items = [
        AnnotationItem(0, 10),
        AnnotationItem(5, 15),
    ]

    positions, densities = compute_feature_density(items, 0, 20, 4)

    assert positions == pytest.approx([2.5, 7.5, 12.5, 17.5])
    assert densities == [1, 2, 1, 0]


@pytest.mark.parametrize(
    "span, expected",
    [(1, 1), (140, 100), (10_000, 5_000), (80_000, 50_000)],
)
def test_scale_ruler_uses_a_readable_125_length(span, expected):
    assert nice_scale_length(span) == expected


@pytest.mark.parametrize(
    "length, expected",
    [(100, "100 bp"), (10_000, "10 kb"), (2_000_000, "2 Mb")],
)
def test_scale_ruler_formats_genomic_units(length, expected):
    assert format_scale_length(length) == expected


def test_scale_ruler_draws_end_caps_and_assembly_label():
    renderer = AlignmentRenderer()
    fig = plt.figure(figsize=(10, 2))

    renderer.draw_scale_bar(fig, 2, 0.10, 0.90, 140, "hg19")

    assert len(fig.lines) == 3
    horizontal = fig.lines[0]
    assert horizontal.get_xdata() == pytest.approx([0.10, 0.10 + 0.80 * 100 / 140])
    assert horizontal.get_ydata()[0] == pytest.approx(horizontal.get_ydata()[1])
    assert {text.get_text() for text in fig.texts} == {"100 bp", "hg19"}
    plt.close(fig)


def test_background_grid_can_be_disabled():
    renderer = AlignmentRenderer(grid_mode="none")
    fig, ax = plt.subplots()

    renderer.draw_background_grid(ax, [20, 40, 60, 80], 0, 100)

    assert not ax.lines
    assert not ax.patches
    plt.close(fig)


def test_major_minor_background_grid_subdivides_coordinate_intervals():
    config = load_config()
    config["styles"]["grid_minor_divisions"] = 2
    renderer = AlignmentRenderer(grid_mode="major_minor", visual_config=config)
    fig, ax = plt.subplots()

    renderer.draw_background_grid(ax, [20, 40, 60, 80], 0, 100)

    major_lines = [line for line in ax.lines if line.get_linestyle() == "-"]
    minor_lines = [line for line in ax.lines if line.get_linestyle() == ":"]
    assert [line.get_xdata()[0] for line in major_lines] == [20, 40, 60, 80]
    assert [line.get_xdata()[0] for line in minor_lines] == [10, 30, 50, 70, 90]
    plt.close(fig)


def test_banded_background_alternates_between_major_boundaries():
    renderer = AlignmentRenderer(grid_mode="bands")
    fig, ax = plt.subplots()

    renderer.draw_background_grid(ax, [20, 40, 60, 80], 0, 100)

    assert len(ax.patches) == 3
    assert len(ax.lines) == 4
    assert [patch.get_x() for patch in ax.patches] == [0, 40, 80]
    assert [patch.get_width() for patch in ax.patches] == [20, 20, 20]
    plt.close(fig)


def test_unknown_background_grid_mode_is_rejected():
    with pytest.raises(ValueError, match="Unknown grid mode"):
        AlignmentRenderer(grid_mode="graph-paper")


def test_highlight_region_is_clipped_and_accepts_chr_aliases():
    renderer = AlignmentRenderer(
        highlight_regions=[
            HighlightRegion("1", 90, 130),
            HighlightRegion("chr2", 110, 150),
        ],
        highlight_color="#12abef",
        highlight_alpha=0.35,
    )
    fig, ax = plt.subplots()

    renderer.draw_highlights(ax, "chr1", 100, 200)

    assert len(ax.patches) == 1
    assert ax.patches[0].get_x() == pytest.approx(100)
    assert ax.patches[0].get_width() == pytest.approx(30)
    assert ax.patches[0].get_alpha() == pytest.approx(0.35)
    assert to_hex(ax.patches[0].get_facecolor()) == "#12abef"
    plt.close(fig)


@pytest.mark.parametrize(
    "kwargs, message",
    [
        ({"highlight_color": "not-a-colour"}, "Invalid highlight color"),
        ({"highlight_alpha": 0}, "Highlight alpha"),
        ({"highlight_alpha": 1.1}, "Highlight alpha"),
        ({"highlight_regions": [HighlightRegion("chr1", 20, 10)]},
         "Highlight region end"),
    ],
)
def test_invalid_highlight_settings_are_rejected(kwargs, message):
    with pytest.raises(ValueError, match=message):
        AlignmentRenderer(**kwargs)


@pytest.mark.parametrize(
    "alignment, expected",
    [
        ("left", (0.01, "left")),
        ("center", (0.50, "center")),
        ("right", (0.99, "right")),
    ],
)
def test_figure_title_position_matches_requested_alignment(alignment, expected):
    renderer = AlignmentRenderer(title_align=alignment)

    assert renderer.figure_title_position() == expected


def test_unknown_figure_title_alignment_is_rejected():
    with pytest.raises(ValueError, match="Unknown title alignment"):
        AlignmentRenderer(title_align="floating")


def test_render_aligns_region_title_and_subtitle_as_one_block(
    tmp_path, monkeypatch,
):
    renderer = AlignmentRenderer(
        title_align="right", show_coverage=False, show_ideogram=False,
        show_legend=False, fig_width=4, dpi=40,
    )
    original_close = plt.close
    monkeypatch.setattr(plt, "close", lambda *_args, **_kwargs: None)
    output = tmp_path / "right-title.png"

    renderer.render(
        rows=[], chrom="chr1", window_start=100, window_end=200,
        reference=None, out_path=str(output), title="Right-aligned subtitle",
    )

    fig = plt.gcf()
    title = next(text for text in fig.texts if text.get_text().startswith("chr1:"))
    subtitle = next(
        text for text in fig.texts if text.get_text() == "Right-aligned subtitle"
    )
    assert title.get_position()[0] == pytest.approx(0.99)
    assert subtitle.get_position()[0] == pytest.approx(0.99)
    assert title.get_ha() == "right"
    assert subtitle.get_ha() == "right"
    original_close(fig)


def test_highlight_is_applied_to_every_data_track_but_not_ideogram(
    tmp_path, monkeypatch,
):
    item = AnnotationItem(100, 110)
    track = LoadedAnnotationTrack(
        "Regions", "bed", "#000000", [item], [[item]], "collapse"
    )
    renderer = AlignmentRenderer(
        show_ideogram=True,
        show_coverage=True,
        show_alignments=True,
        highlight_regions=[HighlightRegion("chr1", 100, 110)],
        fig_width=4,
        dpi=40,
    )
    calls = []
    original = renderer.draw_highlights

    def record_highlight(ax, chrom, start, end):
        calls.append((ax, chrom, start, end))
        original(ax, chrom, start, end)

    monkeypatch.setattr(renderer, "draw_highlights", record_highlight)
    output = tmp_path / "highlight-all-tracks.png"

    renderer.render(
        rows=[], chrom="chr1", window_start=90, window_end=120,
        reference=None, out_path=str(output), genomic_tracks=[track],
        contig_length=1_000, all_reads_for_coverage=[],
    )

    assert output.is_file()
    assert len(calls) == 3  # annotation, coverage, and alignments
    assert all(call[1:] == ("chr1", 90, 120) for call in calls)


def test_peak_track_draws_signal_blocks_and_narrowpeak_summit():
    item = AnnotationItem(100, 180, "peak", value=12.5, summit=130)
    track = LoadedAnnotationTrack(
        "H3K27ac", "narrowpeak", "#7b3294", [item], [[item]], "collapse"
    )
    renderer = AlignmentRenderer()
    fig, ax = plt.subplots()

    renderer.draw_annotation_track(ax, track, 90, 200)

    assert len(ax.patches) == 1
    assert ax.patches[0].get_height() == pytest.approx(12.5)
    assert any(130 in line.get_xdata() for line in ax.lines)
    assert ax.get_ylabel() == "signalValue"
    plt.close(fig)


def test_signal_track_draws_one_continuous_filled_profile_on_fixed_scale():
    items = [
        AnnotationItem(100, 110, value=2.0),
        AnnotationItem(110, 120, value=7.0),
        AnnotationItem(120, 130, value=3.0),
    ]
    track = LoadedAnnotationTrack(
        "CTCF", "signal", "#00695c", items, [items], "collapse"
    )
    config = load_config()
    config["styles"]["signal_y_max"] = 10.0
    renderer = AlignmentRenderer(visual_config=config)
    fig, ax = plt.subplots()

    renderer.draw_annotation_track(ax, track, 95, 135)

    assert not ax.patches
    assert len(ax.collections) == 1
    assert len(ax.lines) == 2
    assert ax.get_ylim()[1] == pytest.approx(11.2)
    assert ax.get_ylabel() == "normalized signal"
    plt.close(fig)


def test_track_only_render_omits_alignment_rows_and_legend(tmp_path, monkeypatch):
    items = [AnnotationItem(100, 110, value=4.0)]
    track = LoadedAnnotationTrack(
        "CTCF", "signal", "#00695c", items, [items], "collapse"
    )
    renderer = AlignmentRenderer(
        show_alignments=False, show_coverage=False, show_ideogram=False,
        fig_width=4, dpi=40,
    )

    def reject_legend(*args, **kwargs):
        raise AssertionError("The alignment legend must not be drawn in track-only mode")

    monkeypatch.setattr(renderer, "draw_legends", reject_legend)
    output = tmp_path / "track-only.png"
    renderer.render(
        rows=[], chrom="chr1", window_start=90, window_end=120,
        reference=None, out_path=str(output), genomic_tracks=[track],
    )

    assert output.is_file()
    assert renderer.legend_height_in == 0


def test_no_legend_reclaims_legend_space_with_alignments_enabled(tmp_path, monkeypatch):
    renderer = AlignmentRenderer(
        show_legend=False, show_coverage=False, show_ideogram=False,
        fig_width=4, dpi=40,
    )

    def reject_legend(*args, **kwargs):
        raise AssertionError("The legend must not be drawn with show_legend=False")

    monkeypatch.setattr(renderer, "draw_legends", reject_legend)
    output = tmp_path / "no-legend.png"
    renderer.render(
        rows=[], chrom="chr1", window_start=90, window_end=120,
        reference=None, out_path=str(output),
    )

    assert output.is_file()
    assert renderer.show_alignments
    assert not renderer.show_legend
    assert renderer.legend_height_in == 0
    assert renderer.legend_margin_in == pytest.approx(
        renderer.legend_tick_clearance_in
    )


def test_only_bottom_coordinate_axis_has_x_tick_marks(tmp_path, monkeypatch):
    renderer = AlignmentRenderer(
        show_alignments=False, show_coverage=True, show_ideogram=True,
        fig_width=4, dpi=40,
    )
    original_close = plt.close
    monkeypatch.setattr(plt, "close", lambda *_args, **_kwargs: None)
    output = tmp_path / "axis-ticks.png"

    renderer.render(
        rows=[], chrom="chr1", window_start=100, window_end=200,
        reference=None, out_path=str(output), contig_length=1_000,
        all_reads_for_coverage=[],
    )

    fig = plt.gcf()
    ideogram_ax, coordinate_ax = fig.axes
    assert not any(
        tick.tick1line.get_visible() or tick.tick2line.get_visible()
        for tick in ideogram_ax.xaxis.get_major_ticks()
    )
    assert any(
        tick.tick1line.get_visible()
        for tick in coordinate_ax.xaxis.get_major_ticks()
    )
    original_close(fig)


def test_density_track_draws_compact_filled_histogram():
    items = [AnnotationItem(100, 150), AnnotationItem(125, 175)]
    track = LoadedAnnotationTrack(
        "DNase", "narrowpeak", "#2c7fb8", items, [items], "density"
    )
    renderer = AlignmentRenderer()
    fig, ax = plt.subplots()

    renderer.draw_annotation_track(ax, track, 90, 190)

    assert len(ax.collections) == 1
    assert ax.get_ylabel() == "features/bin"
    assert any(text.get_text() == "density" for text in ax.texts)
    plt.close(fig)


def test_annotation_height_override_wins_over_type_default():
    item = AnnotationItem(100, 150)
    track = LoadedAnnotationTrack(
        "Custom", "bed", "#000000", [item], [[item]], height_in=0.72
    )
    renderer = AlignmentRenderer()

    assert renderer.annotation_track_height(track) == pytest.approx(0.72)
    assert renderer.styles["coverage_track_height_in"] == pytest.approx(1.40)


def test_compute_coverage_counts_match_blocks_but_not_deletions():
    reads = [
        SimpleNamespace(blocks=[
            CigarBlock("M", 100, 0, 3),
            CigarBlock("D", 103, 3, 2),
            CigarBlock("M", 105, 3, 2),
        ]),
        SimpleNamespace(blocks=[CigarBlock("M", 101, 0, 5)]),
    ]
    assert compute_coverage(reads, 100, 107) == [1, 2, 2, 1, 1, 2, 1]


def test_binned_coverage_returns_mean_depth_without_per_base_output():
    reads = [
        SimpleNamespace(blocks=[CigarBlock("M", 0, 0, 4)]),
        SimpleNamespace(blocks=[CigarBlock("M", 2, 0, 6)]),
    ]

    edges, depth = compute_binned_coverage(reads, 0, 10, 2)

    assert edges == pytest.approx([0, 5, 10])
    assert depth == pytest.approx([1.4, 0.6])


def test_coverage_colors_only_snvs_above_vaf_threshold():
    reads = []
    for index in range(5):
        mismatches = []
        if index < 2:
            mismatches.append((100, "A"))
        if index == 2:
            mismatches.append((100, "C"))
        reads.append(SimpleNamespace(
            blocks=[CigarBlock("M", 100, 0, 1)], mismatches=mismatches,
            mismatch_details=[(position, base, 35) for position, base in mismatches],
            query_sequence=mismatches[0][1] if mismatches else "G",
            query_qualities=[35], mapq=60, is_reverse=index % 2 == 1,
        ))

    assert compute_snv_counts(reads, 100, 101) == {100: {"A": 2, "C": 1}}

    renderer = AlignmentRenderer(
        coverage_vaf_threshold=0.20, show_variant_counts=True
    )
    fig, ax = plt.subplots()
    renderer.draw_coverage_track(ax, reads, 100, 101)

    assert len(ax.patches) == 2  # grey depth plus the 40% A allele; 20% C is hidden
    assert ax.patches[1].get_height() == 2
    assert to_hex(ax.patches[1].get_facecolor()) == DEFAULT_BASE_COLORS["A"]
    assert any(
        text.get_text() == "A 2/5 40% F1/R1 BQ35 MQ60"
        for text in ax.texts
    )
    plt.close(fig)


def test_snv_evidence_applies_quality_filters_and_tracks_strand_means():
    reads = [
        SimpleNamespace(
            blocks=[CigarBlock("M", 100, 0, 1)], query_sequence="G",
            query_qualities=[35], mapq=60, is_reverse=False,
            mismatch_details=[(100, "G", 35)],
        ),
        SimpleNamespace(
            blocks=[CigarBlock("M", 100, 0, 1)], query_sequence="G",
            query_qualities=[25], mapq=40, is_reverse=True,
            mismatch_details=[(100, "G", 25)],
        ),
        SimpleNamespace(
            blocks=[CigarBlock("M", 100, 0, 1)], query_sequence="G",
            query_qualities=[10], mapq=60, is_reverse=False,
            mismatch_details=[(100, "G", 10)],
        ),
    ]

    depth, evidence = compute_snv_evidence(
        reads, 100, 101, min_baseq=20, min_mapq=30
    )
    allele = evidence[100]["G"]

    assert depth == [2]
    assert (allele.count, allele.forward, allele.reverse) == (2, 1, 1)
    assert allele.base_quality_sum / allele.count == 30
    assert allele.mapq_sum / allele.count == 50

    sparse_depth, sparse_evidence = compute_sparse_snv_evidence(
        reads, 100, 101, min_baseq=20, min_mapq=30
    )
    sparse_allele = sparse_evidence[100]["G"]
    assert sparse_depth == {100: 2}
    assert (sparse_allele.count, sparse_allele.forward, sparse_allele.reverse) == (2, 1, 1)


def test_wide_coverage_uses_binned_single_artist():
    read = SimpleNamespace(
        blocks=[CigarBlock("M", 100, 0, 150)],
        mismatches=[], mismatch_details=[], query_sequence="A" * 150,
        query_qualities=[35] * 150, mapq=60, is_reverse=False,
    )
    renderer = AlignmentRenderer()
    fig, ax = plt.subplots(figsize=(6, 1), dpi=100)

    renderer.draw_coverage_track(ax, [read], 0, 100_000)

    assert len(ax.patches) == 1
    assert len(ax.patches[0].get_path().vertices) < 2_000
    assert any(text.get_text().endswith("bp/bin (mean)") for text in ax.texts)
    plt.close(fig)


def test_squish_rows_are_shorter_than_expanded_rows():
    expanded = AlignmentRenderer(display_mode="expand")
    squished = AlignmentRenderer(display_mode="squish")
    assert squished.row_height_in < expanded.row_height_in
    assert squished.row_height_in <= expanded.row_height_in / 4
    assert squished.row_margin < expanded.row_margin
    assert squished.styles["squish_alignment_edge_width"] == 0


def test_close_zoom_labels_each_softclipped_nucleotide():
    renderer = AlignmentRenderer(display_mode="expand", shade_by_mapq=False)
    read = SimpleNamespace(
        ref_start=100, ref_end=104, pair_category="normal", mate_chrom="chr1",
        is_secondary=False, is_duplicate=False, mapq=60,
        blocks=[
            CigarBlock("S", 100, 0, 4),
            CigarBlock("M", 100, 4, 4),
            CigarBlock("S", 104, 8, 4),
        ],
        mismatches=[], query_sequence="ACGTAAAATGCA",
    )
    fig, ax = plt.subplots(figsize=(4, 1), dpi=100)
    ax.set_xlim(90, 110)

    renderer.draw_read(ax, read, y0=0.1, h=0.8, render_base_detail=True)

    assert [text.get_text() for text in ax.texts] == list("ACGTTGCA")
    assert [text.get_position()[0] for text in ax.texts] == pytest.approx(
        [96.5, 97.5, 98.5, 99.5, 104.5, 105.5, 106.5, 107.5]
    )
    assert [to_hex(patch.get_facecolor()) for patch in ax.patches[:4]] == [
        renderer.alignment_colors["normal"]
    ] * 4
    assert [to_hex(text.get_color()) for text in ax.texts[:4]] == [
        DEFAULT_BASE_COLORS[base] for base in "ACGT"
    ]
    plt.close(fig)


def test_wider_view_keeps_softclip_cells_but_omits_letters():
    renderer = AlignmentRenderer(display_mode="expand", shade_by_mapq=False)
    read = SimpleNamespace(
        ref_start=100, ref_end=104, pair_category="normal", mate_chrom="chr1",
        is_secondary=False, is_duplicate=False, mapq=60,
        blocks=[CigarBlock("S", 100, 0, 4), CigarBlock("M", 100, 4, 4)],
        mismatches=[], query_sequence="ACGTAAAA",
    )
    fig, ax = plt.subplots(figsize=(4, 1), dpi=100)
    ax.set_xlim(0, 200)

    renderer.draw_read(ax, read, y0=0.1, h=0.8, render_base_detail=True)

    assert len(ax.patches) == 5
    assert not ax.texts
    plt.close(fig)


def test_close_zoom_draws_igv_style_insertion_marker():
    renderer = AlignmentRenderer(display_mode="expand", shade_by_mapq=False)
    read = SimpleNamespace(
        ref_start=100, ref_end=120, pair_category="normal", mate_chrom="chr1",
        is_secondary=False, is_duplicate=False, mapq=60,
        blocks=[
            CigarBlock("M", 100, 0, 10),
            CigarBlock("I", 110, 10, 4),
            CigarBlock("M", 110, 14, 10),
        ],
        mismatches=[], query_sequence="A" * 24,
    )
    fig, ax = plt.subplots(figsize=(8, 1), dpi=100)
    ax.set_xlim(90, 130)

    renderer.draw_read(ax, read, y0=0.1, h=0.8, render_base_detail=True)

    marker = ax.patches[1]
    assert marker.get_x() + marker.get_width() / 2 == pytest.approx(110)
    assert marker.get_width() == pytest.approx(0.4)
    assert to_hex(marker.get_facecolor()) == DEFAULT_VISUAL_COLORS["insertion"]
    assert [text.get_text() for text in ax.texts] == ["I"]
    assert to_hex(ax.texts[0].get_color()) == DEFAULT_VISUAL_COLORS["contrast_edge"]
    plt.close(fig)


def test_squish_reads_hide_outlines_and_per_read_event_labels():
    renderer = AlignmentRenderer(display_mode="squish", shade_by_mapq=False)
    read = SimpleNamespace(
        ref_start=100, ref_end=125, pair_category="normal", mate_chrom="chr1",
        is_secondary=False, is_duplicate=False, mapq=60,
        blocks=[
            CigarBlock("M", 100, 0, 10),
            CigarBlock("D", 110, 10, 5),
            CigarBlock("I", 115, 10, 4),
            CigarBlock("M", 115, 14, 10),
        ],
        mismatches=[], query_sequence="A" * 24,
    )
    fig, ax = plt.subplots()

    renderer.draw_read(ax, read, y0=0.08, h=0.84, render_base_detail=False)

    assert all(patch.get_linewidth() == 0 for patch in ax.patches)
    assert not ax.texts
    plt.close(fig)


def test_expanded_reads_have_no_outline_by_default():
    renderer = AlignmentRenderer(display_mode="expand", shade_by_mapq=False)
    read = SimpleNamespace(
        ref_start=100, ref_end=120, pair_category="normal", mate_chrom="chr1",
        is_secondary=False, is_duplicate=False, mapq=60,
        blocks=[CigarBlock("M", 100, 0, 20)],
        mismatches=[], query_sequence="A" * 20,
    )
    fig, ax = plt.subplots()

    renderer.draw_read(ax, read, y0=0.1, h=0.8, render_base_detail=False)

    assert len(ax.patches) == 1
    assert ax.patches[0].get_linewidth() == 0
    assert to_hex(ax.patches[0].get_facecolor()) == renderer.alignment_colors["normal"]
    plt.close(fig)


def test_indel_length_labels_are_opt_in():
    read = SimpleNamespace(
        ref_start=100, ref_end=125, pair_category="normal", mate_chrom="chr1",
        is_secondary=False, is_duplicate=False, mapq=60,
        blocks=[
            CigarBlock("M", 100, 0, 10),
            CigarBlock("D", 110, 10, 5),
            CigarBlock("I", 115, 10, 4),
            CigarBlock("M", 115, 14, 10),
        ],
        mismatches=[], query_sequence="A" * 24,
    )
    fig, (default_ax, labelled_ax) = plt.subplots(nrows=2)

    AlignmentRenderer(shade_by_mapq=False).draw_read(
        default_ax, read, y0=0.1, h=0.8, render_base_detail=False
    )
    AlignmentRenderer(
        shade_by_mapq=False, show_indel_lengths=True
    ).draw_read(
        labelled_ax, read, y0=0.1, h=0.8, render_base_detail=False
    )

    assert [text.get_text() for text in default_ax.texts] == ["I"]
    assert [text.get_text() for text in labelled_ax.texts] == ["5", "I", "+4"]
    plt.close(fig)


def test_split_read_junction_is_a_thin_solid_line():
    renderer = AlignmentRenderer(shade_by_mapq=False)
    read = SimpleNamespace(
        ref_start=100, ref_end=140, pair_category="normal", mate_chrom="chr1",
        is_secondary=False, is_duplicate=False, mapq=60,
        blocks=[
            CigarBlock("M", 100, 0, 10),
            CigarBlock("N", 110, 10, 20),
            CigarBlock("M", 130, 10, 10),
        ],
        mismatches=[], query_sequence="A" * 20,
    )
    fig, ax = plt.subplots()

    renderer.draw_read(ax, read, y0=0.1, h=0.8, render_base_detail=False)

    assert len(ax.lines) == 1
    assert ax.lines[0].get_linestyle() == "-"
    assert ax.lines[0].get_linewidth() == pytest.approx(0.55)
    plt.close(fig)


def test_deletion_connector_uses_thin_default_line_width():
    renderer = AlignmentRenderer(shade_by_mapq=False)
    read = SimpleNamespace(
        ref_start=100, ref_end=1120, pair_category="normal", mate_chrom="chr1",
        is_secondary=False, is_duplicate=False, mapq=60,
        blocks=[
            CigarBlock("M", 100, 0, 10),
            CigarBlock("D", 110, 10, 1000),
            CigarBlock("M", 1110, 10, 10),
        ],
        mismatches=[], query_sequence="A" * 20,
    )
    fig, ax = plt.subplots()

    renderer.draw_read(ax, read, y0=0.1, h=0.8, render_base_detail=False)

    assert len(ax.lines) == 1
    assert ax.lines[0].get_linewidth() == pytest.approx(0.65)
    plt.close(fig)


def test_long_labels_are_ellipsized_to_the_available_lane():
    assert ellipsize("Candidate regions", 10) == "Candidate…"
    assert ellipsize("short", 10) == "short"


@pytest.mark.parametrize(
    "ticks, window_length, expected_labels, expected_unit",
    [
        (
            [101_867_480, 101_867_500, 101_867_520], 140,
            ["101,867,480", "101,867,500", "101,867,520"], "bp",
        ),
        (
            [101_860_000, 101_865_000, 101_870_000], 20_000,
            ["101,860", "101,865", "101,870"], "kb",
        ),
        (
            [101_000_000, 101_200_000, 101_400_000], 1_400_000,
            ["101.0", "101.2", "101.4"], "Mb",
        ),
        (
            [3_000_000_000, 3_100_000_000], 1_100_000_000,
            ["3.0", "3.1"], "Gb",
        ),
    ],
)
def test_genomic_ticks_use_explicit_units_not_exponents(
    ticks, window_length, expected_labels, expected_unit,
):
    labels, unit = genomic_tick_labels(ticks, window_length)

    assert labels == expected_labels
    assert unit == expected_unit
    assert not any("e" in label.lower() for label in labels)


def test_legend_clusters_related_terms_by_topic():
    renderer = AlignmentRenderer(fig_width=14)
    fig = plt.figure(figsize=(14, 2))
    legends = renderer.draw_legends(fig, fig_height=2)

    assert [legend.get_title().get_text() for legend in legends] == [
        "Alignment", "Read events", "Insert size", "Base identity",
    ]
    assert [text.get_text() for text in legends[0].get_texts()] == [
        "Normal / concordant", "FF (same strand)", "RR (same strand)",
        "Reverted (RF)", "Inter-chromosomal (mate colour)",
    ]
    assert [text.get_text() for text in legends[1].get_texts()] == [
        "Insertion", "Deletion",
    ]
    assert [text.get_text() for text in legends[3].get_texts()] == list("ACGT")
    legend_ax = fig.axes[-1]
    assert len(legend_ax.patches) == 4
    assert len(legend_ax.lines) == 0
    bounds = []
    for patch in legend_ax.patches:
        bounds.append((patch.get_x(), patch.get_x() + patch.get_width()))
    for left, right in zip(bounds, bounds[1:]):
        assert left[1] < right[0]
    plt.close(fig)


def test_alignment_legend_lists_at_most_two_mate_chromosomes_explicitly():
    renderer = AlignmentRenderer(fig_width=14)
    renderer.interchrom_mate_colors = {
        "chr2": "#ce3d32",
        "chr1": "#5050ff",
    }
    fig = plt.figure(figsize=(14, 2))

    legends = renderer.draw_legends(fig, fig_height=2)

    labels = [text.get_text() for text in legends[0].get_texts()]
    assert labels[-2:] == ["Mate chr1", "Mate chr2"]
    plt.close(fig)


def test_alignment_legend_collapses_many_mate_chromosomes_to_one_entry():
    renderer = AlignmentRenderer(fig_width=14)
    renderer.interchrom_mate_colors = {
        f"chr{chromosome}": f"#{chromosome:02x}55aa"
        for chromosome in range(1, 9)
    }
    fig = plt.figure(figsize=(14, 2))

    legends = renderer.draw_legends(fig, fig_height=2)

    labels = [text.get_text() for text in legends[0].get_texts()]
    assert labels == [
        "Normal / concordant", "FF (same strand)", "RR (same strand)",
        "Reverted (RF)", "Inter-chromosomal (8 chromosomes)",
    ]
    assert not any(label.startswith("Mate chr") for label in labels)
    plt.close(fig)


def test_haplotype_view_colours_reads_and_replaces_pair_legend_compartment():
    renderer = AlignmentRenderer(haplotype_view="color")
    read = SimpleNamespace(
        haplotype="2", pair_category="large_insert", mate_chrom="chr1",
        is_secondary=False, is_duplicate=False, mapq=60,
    )
    color, alpha = renderer.read_style(read)
    assert color == haplotype_color("2")
    assert alpha == pytest.approx(0.9)

    fig = plt.figure(figsize=(14, 2))
    legends = renderer.draw_legends(fig, fig_height=2)
    assert [legend.get_title().get_text() for legend in legends] == [
        "Alignment", "Read events", "Haplotype", "Base identity",
    ]
    assert [text.get_text() for text in legends[2].get_texts()] == [
        "HP 1", "HP 2", "Other HP", "Untagged",
    ]
    plt.close(fig)


def test_split_haplotype_lanes_show_hp_and_phase_set_labels():
    renderer = AlignmentRenderer(haplotype_view="split")
    rows = [
        [SimpleNamespace(haplotype="1", phase_set="100")],
        [SimpleNamespace(haplotype="1", phase_set="100")],
        [SimpleNamespace(haplotype="2", phase_set="200")],
        [SimpleNamespace(haplotype=None, phase_set=None)],
    ]
    fig, ax = plt.subplots(figsize=(8, 2))
    ax.set_ylim(4, 0)

    renderer.draw_haplotype_lanes(ax, rows)

    assert [text.get_text() for text in ax.texts] == [
        "HP 1 · PS 100", "HP 2 · PS 200", "Untagged",
    ]
    assert len(ax.lines) == 2
    assert len(ax.patches) == 2
    plt.close(fig)


def test_generic_tag_view_colours_reads_and_builds_dynamic_legend():
    renderer = AlignmentRenderer(
        read_tag="RG", tag_view="color", tag_label="Library",
        tag_colors={"tumour": "#445566"},
    )
    read = SimpleNamespace(
        tag_value="tumour", pair_category="large_insert", mate_chrom="chr1",
        is_secondary=False, is_duplicate=False, mapq=60,
    )

    color, alpha = renderer.read_style(read)

    assert color == "#445566"
    assert tag_color("tumour", {"tumour": "#445566"}) == "#445566"
    assert alpha == pytest.approx(0.9)
    fig = plt.figure(figsize=(14, 2))
    legends = renderer.draw_legends(fig, fig_height=2)
    assert [legend.get_title().get_text() for legend in legends] == [
        "Alignment", "Read events", "Library", "Base identity",
    ]
    assert [text.get_text() for text in legends[2].get_texts()] == ["tumour (n=1)"]
    plt.close(fig)


def test_split_generic_tag_lanes_use_tag_label_and_untagged_value():
    renderer = AlignmentRenderer(
        read_tag="RG", tag_view="split", tag_label="Library"
    )
    rows = [
        [SimpleNamespace(tag_value="A")],
        [SimpleNamespace(tag_value="A")],
        [SimpleNamespace(tag_value="B")],
        [SimpleNamespace(tag_value=None)],
    ]
    fig, ax = plt.subplots(figsize=(8, 2))
    ax.set_ylim(4, 0)

    renderer.draw_haplotype_lanes(ax, rows)

    assert [text.get_text() for text in ax.texts] == [
        "Library=A (n=2)", "Library=B (n=1)", "Library=untagged (n=1)",
    ]
    assert len(ax.lines) == 2
    assert len(ax.patches) == 2
    plt.close(fig)


def test_high_cardinality_tag_legend_is_bounded():
    renderer = AlignmentRenderer(read_tag="CB", tag_view="color")
    for index in range(12):
        renderer.read_style(SimpleNamespace(
            tag_value=f"cell-{index:02d}", pair_category="normal", mate_chrom="chr1",
            is_secondary=False, is_duplicate=False, mapq=60,
        ))
    fig = plt.figure(figsize=(14, 2))

    legends = renderer.draw_legends(fig, fig_height=2)
    labels = [text.get_text() for text in legends[2].get_texts()]

    assert len(labels) == 8
    assert "5 more values" in labels
    plt.close(fig)


@pytest.mark.parametrize("fig_width", [5, 8, 14])
def test_legend_keeps_rendered_text_clear_of_plot_and_coordinate_labels(fig_width):
    renderer = AlignmentRenderer(fig_width=fig_width, view_as_pairs=True)
    fig_height = 6
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    fig.subplots_adjust(
        left=0.08, right=0.92, top=0.95,
        bottom=renderer.legend_margin_in / fig_height,
    )
    legends = renderer.draw_legends(fig, fig_height, 0.08, 0.92)
    renderer.separate_legend_from_plots(fig, [ax], legends)

    canvas_renderer = fig.canvas.get_renderer()
    legend_top = max(
        legends[0].axes.get_window_extent(canvas_renderer).y1,
        max(legend.get_window_extent(canvas_renderer).y1 for legend in legends),
    )
    plot_bottom = ax.get_tightbbox(canvas_renderer).y0
    assert plot_bottom >= legend_top + renderer.legend_plot_gap_in * fig.dpi - 1
    plt.close(fig)


def test_small_window_reference_track_draws_coloured_base_cells_and_letters():
    renderer = AlignmentRenderer(max_reference_span=4)
    sequence = "ACGT"
    reference = SimpleNamespace(
        available=True,
        base_at=lambda position: sequence[position],
    )
    fig, ax = plt.subplots(figsize=(4, 1))
    ax.set_xlim(0, 4)

    renderer.draw_reference_track(ax, reference, 0, 4, available_width_in=4)

    assert len(ax.patches) == 4
    assert [text.get_text() for text in ax.texts] == ["A", "C", "G", "T"]
    plt.close(fig)


def test_view_as_pairs_draws_a_link_between_visible_primary_mates():
    renderer = AlignmentRenderer(view_as_pairs=True, shade_by_mapq=False)
    common = dict(
        query_name="pair", reference_name="chr1", mate_chrom="chr1",
        is_paired=True, is_secondary=False, is_supplementary=False,
        mate_is_unmapped=False, is_duplicate=False, pair_category="normal",
        mapq=60, blocks=[], gap_label=lambda: "",
    )
    left = SimpleNamespace(ref_start=10, ref_end=20, **common)
    right = SimpleNamespace(ref_start=40, ref_end=50, **common)
    fig, ax = plt.subplots()

    renderer.draw_alignment_row(
        ax, [left, right], y0=0.1, h=0.8,
        render_base_detail=False, layout="expand",
    )

    assert len(ax.lines) == 1
    assert list(ax.lines[0].get_xdata()) == [20, 40]
    assert len(ax.patches) == 2
    plt.close(fig)


def test_base_sort_highlights_the_alternative_allele_cell():
    renderer = AlignmentRenderer(
        sort_base_position=10, sort_reference_base="C", shade_by_mapq=False
    )
    read = SimpleNamespace(
        ref_start=5, ref_end=15, pair_category="normal", mate_chrom="chr1",
        is_secondary=False, is_duplicate=False, mapq=60,
        blocks=[CigarBlock("M", 5, 0, 10)], mismatches=[],
        query_sequence="CCCCCACCCC",
        base_at=lambda position: "A" if position == 10 else "C",
    )
    fig, ax = plt.subplots()

    renderer.draw_read(ax, read, y0=0.1, h=0.8, render_base_detail=True)

    assert len(ax.patches) == 2
    assert ax.patches[-1].get_x() == 10
    assert to_hex(ax.patches[-1].get_facecolor()) == DEFAULT_BASE_COLORS["A"]
    plt.close(fig)


def test_center_guide_is_hidden_by_default_and_centered_when_enabled():
    fig, (hidden_ax, visible_ax) = plt.subplots(nrows=2)
    AlignmentRenderer().draw_center_guide(hidden_ax, 100, 200)
    AlignmentRenderer(show_center_guide=True).draw_center_guide(
        visible_ax, 100, 200
    )

    assert len(hidden_ax.lines) == 0
    assert len(visible_ax.lines) == 1
    assert list(visible_ax.lines[0].get_xdata()) == [150, 150]
    assert visible_ax.lines[0].get_linestyle() == "--"
    plt.close(fig)


def test_sashimi_counts_cigar_skips_and_can_split_strands():
    reads = [
        SimpleNamespace(strand="+", deletions=[(120, 80, True)]),
        SimpleNamespace(strand="+", deletions=[(120, 80, True)]),
        SimpleNamespace(strand="-", deletions=[(120, 80, True)]),
        SimpleNamespace(strand="+", deletions=[(140, 10, False)]),
    ]

    assert compute_splice_junctions(reads) == {(120, 200, "."): 3}
    assert compute_splice_junctions(reads, "split") == {
        (120, 200, "+"): 2,
        (120, 200, "-"): 1,
    }


def test_sashimi_draws_supported_count_labelled_arcs():
    reads = [
        SimpleNamespace(strand="+", deletions=[(120, 80, True)]),
        SimpleNamespace(strand="+", deletions=[(120, 80, True)]),
        SimpleNamespace(strand="-", deletions=[(120, 80, True)]),
    ]
    renderer = AlignmentRenderer(
        show_sashimi=True, min_junction_reads=1, sashimi_strand="split"
    )
    fig, ax = plt.subplots()

    renderer.draw_sashimi_track(ax, reads, 100, 220)

    assert len(ax.patches) == 2
    labels = set()
    for text in ax.texts:
        labels.add(text.get_text())
    assert {"splice junctions", "+2", "-1"}.issubset(labels)
    assert ax.get_ylim() == pytest.approx((-1.05, 1.05))
    plt.close(fig)


def test_ideogram_marks_the_window_in_red():
    renderer = AlignmentRenderer()
    fig, ax = plt.subplots()
    renderer.draw_ideogram(ax, "chr1", 100, 200, 1_000)
    fig.canvas.draw()

    assert len(ax.patches) == 2
    assert ax.patches[0].get_x() == 100
    assert ax.patches[0].get_width() == 100
    chromosome_bounds = ax.patches[0].get_window_extent(fig.canvas.get_renderer())
    axes_bounds = ax.get_window_extent()
    assert chromosome_bounds.x0 == pytest.approx(axes_bounds.x0)
    assert chromosome_bounds.x1 == pytest.approx(axes_bounds.x1)
    assert to_hex(ax.patches[1].get_facecolor()) == DEFAULT_VISUAL_COLORS["ideogram_window"]
    assert {text.get_text() for text in ax.texts} == {"chr1", "0.0 Mb"}
    plt.close(fig)


def test_ideogram_keeps_a_visible_marker_for_a_tiny_genomic_window():
    renderer = AlignmentRenderer()
    fig, ax = plt.subplots()
    renderer.draw_ideogram(ax, "chr1", 100, 200, 250_000_000)

    chromosome = ax.patches[0]
    marker = ax.patches[-1]
    assert marker.get_width() == pytest.approx(chromosome.get_width() * 0.004)
    assert to_hex(marker.get_facecolor()) == DEFAULT_VISUAL_COLORS["ideogram_window"]
    assert "window" not in {text.get_text() for text in ax.texts}
    plt.close(fig)


def test_ideogram_draws_cytobands_and_centromere():
    renderer = AlignmentRenderer()
    fig, ax = plt.subplots()
    bands = [
        Cytoband("chr1", 0, 400, "p11", "gneg"),
        Cytoband("chr1", 400, 500, "p10", "acen"),
        Cytoband("chr1", 500, 600, "q10", "acen"),
        Cytoband("chr1", 600, 1_000, "q11", "gpos100"),
    ]
    renderer.draw_ideogram(ax, "chr1", 100, 200, 1_000, bands)

    assert len(ax.patches) == 8  # base, four bands, bridge, outline, marker
    assert sum(isinstance(patch, Polygon) for patch in ax.patches) == 2
    chromosome_vertices = ax.patches[0].get_xy()
    neck_y = [
        y for x, y in chromosome_vertices if x == pytest.approx(150)
    ]
    assert len(neck_y) == 2
    assert min(neck_y) > 0.28
    assert max(neck_y) < 0.72
    assert min(neck_y) < 0.5 < max(neck_y)
    assert to_hex(ax.patches[5].get_facecolor()) == "#b84b4b"
    assert to_hex(ax.patches[-1].get_facecolor()) == DEFAULT_VISUAL_COLORS["ideogram_window"]
    plt.close(fig)


@pytest.mark.parametrize("strand, marker", [("+", ">"), ("-", "<")])
def test_every_gene_intron_has_strand_orientation_arrow(strand, marker):
    renderer = AlignmentRenderer()
    item = AnnotationItem(
        90, 210, "TX1", strand,
        blocks=[(100, 120), (180, 190), (191, 200)],
    )
    track = LoadedAnnotationTrack("Genes", "gtf", "#17217a", [item], [[item]])
    fig, ax = plt.subplots(figsize=(8, 1))
    ax.set_xlim(90, 210)

    renderer.draw_annotation_track(ax, track, 90, 210)

    arrow_positions = []
    for line in ax.lines:
        if line.get_marker() == marker:
            arrow_positions.append(line.get_xdata()[0])
    assert any(90 < position < 100 for position in arrow_positions)
    assert any(120 < position < 180 for position in arrow_positions)
    assert any(190 < position < 191 for position in arrow_positions)
    assert any(200 < position < 210 for position in arrow_positions)
    opposite_marker = "<" if marker == ">" else ">"
    assert not any(line.get_marker() == opposite_marker for line in ax.lines)
    plt.close(fig)


def test_gene_orientation_arrows_use_readable_default_spacing_and_size():
    renderer = AlignmentRenderer()
    item = AnnotationItem(
        100, 200, "TX1", "+", blocks=[(100, 110), (190, 200)]
    )
    track = LoadedAnnotationTrack("Genes", "gtf", "#17217a", [item], [[item]])
    fig, ax = plt.subplots(figsize=(8, 1), dpi=100)
    ax.set_xlim(100, 200)

    renderer.draw_annotation_track(ax, track, 100, 200)

    arrows = []
    for line in ax.lines:
        if line.get_marker() == ">":
            arrows.append(line)
    expected_maximum = int(ax.get_window_extent().width / 28)
    assert 1 < len(arrows) <= expected_maximum
    assert all(
        line.get_markersize() == pytest.approx(3.1) for line in arrows
    )
    plt.close(fig)


def test_primary_isoform_label_has_visible_marker():
    renderer = AlignmentRenderer()
    item = AnnotationItem(
        100, 160, "GENE1", "+", blocks=[(100, 160)],
        group="g1", group_label="GENE1", transcript_label="TX1",
        primary_rank=1, primary_label="MANE Select",
    )
    track = LoadedAnnotationTrack(
        "Genes", "gtf", "#17217a", [item], [[item]], display_mode="expand"
    )
    fig, ax = plt.subplots(figsize=(8, 1))
    ax.set_xlim(90, 170)

    renderer.draw_annotation_track(ax, track, 90, 170)

    assert any(text.get_text() == "GENE1 · TX1 ★" for text in ax.texts)
    plt.close(fig)


def test_expanded_gene_label_falls_back_to_gene_id():
    renderer = AlignmentRenderer()
    item = AnnotationItem(
        100, 160, "tx1", "+", blocks=[(100, 160)],
        group="gene_id_1", transcript_label="tx1",
    )
    track = LoadedAnnotationTrack(
        "Genes", "gtf", "#17217a", [item], [[item]], display_mode="expand"
    )
    fig, ax = plt.subplots(figsize=(8, 1))
    ax.set_xlim(90, 170)

    renderer.draw_annotation_track(ax, track, 90, 170)

    assert any(text.get_text() == "gene_id_1 · tx1" for text in ax.texts)
    plt.close(fig)


def test_cnv_track_draws_log2_segments_around_zero_with_gain_loss_colours():
    renderer = AlignmentRenderer()
    loss = AnnotationItem(100, 130, value=-0.6, sample="Tumour")
    gain = AnnotationItem(130, 170, value=0.8, sample="Tumour")
    track = LoadedAnnotationTrack(
        "Tumour CNV", "seg", "#555555", [loss, gain], [[loss, gain]],
        color_by_sign=True,
    )
    fig, ax = plt.subplots(figsize=(8, 2))
    ax.set_xlim(90, 180)

    renderer.draw_annotation_track(ax, track, 90, 180)

    assert len(ax.patches) == 2
    assert [to_hex(patch.get_facecolor()) for patch in ax.patches] == [
        DEFAULT_VISUAL_COLORS["cnv_loss"], DEFAULT_VISUAL_COLORS["cnv_gain"],
    ]
    assert [line.get_ydata()[0] for line in ax.lines] == [0, -0.6, 0.8]
    assert [tick.get_text() for tick in ax.get_yticklabels()] == ["-1", "0", "1"]
    assert any(text.get_text() == "Tumour" for text in ax.texts)
    plt.close(fig)


def test_baf_track_draws_heterozygous_snvs_around_half_baseline():
    renderer = AlignmentRenderer()
    low = AnnotationItem(100, 101, "rs1", value=0.18, sample="Tumour")
    balanced = AnnotationItem(130, 131, "rs2", value=0.52, sample="Tumour")
    high = AnnotationItem(160, 161, "rs3", value=0.84, sample="Tumour")
    track = LoadedAnnotationTrack(
        "Tumour BAF", "baf", "#7a1f5c",
        [low, balanced, high], [[low, balanced, high]],
    )
    fig, ax = plt.subplots(figsize=(8, 2))
    ax.set_xlim(90, 180)

    renderer.draw_annotation_track(ax, track, 90, 180)

    assert len(ax.collections) == 1
    assert ax.collections[0].get_offsets().tolist() == [
        [100.5, 0.18], [130.5, 0.52], [160.5, 0.84],
    ]
    assert list(ax.lines[0].get_ydata()) == [0.5, 0.5]
    assert [tick.get_text() for tick in ax.get_yticklabels()] == ["0.0", "0.5", "1.0"]
    assert ax.get_ylabel() == "BAF"
    plt.close(fig)
