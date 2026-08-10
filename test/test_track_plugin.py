import os
import sys
from types import MappingProxyType

import matplotlib.pyplot as plt
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from locus_snap import (
    PluginTrackSource,
    TrackCanvas,
    TrackPlugin,
    TrackPluginError,
    TrackRegion,
    load_track_plugin,
)
from locus_snap.config import load_config
from locus_snap.render import AlignmentRenderer
from locus_snap.snapshot import BamSnapshot
from locus_snap.track_plugin import build_plugin_track_sources
import locus_snap.track_plugin as plugin_api


class RecordingPlugin(TrackPlugin):
    name = "recording"
    default_label = "External values"
    default_height_in = 0.72
    default_color = "#123456"

    def fetch(self, region, options):
        assert isinstance(options, MappingProxyType)
        return {"x": [region.start, region.end], "y": [0.2, 0.8]}

    def render(self, canvas, payload, region, options):
        assert canvas.region == region
        assert isinstance(options, MappingProxyType)
        canvas.set_y_limits(0, 1)
        canvas.set_y_ticks([0, 1], ["0", "1"])
        canvas.fill_between(payload["x"], payload["y"], alpha=0.2)
        canvas.line(payload["x"], payload["y"])


def test_plugin_source_fetches_read_only_options_and_metadata():
    source = PluginTrackSource(
        RecordingPlugin(), {"sample": "tumour"},
        label="Plugin values", height_in=1.1, color="#abcdef",
    )

    track = source.fetch("chr1", 100, 200)

    assert track.label == "Plugin values"
    assert track.height_in == pytest.approx(1.1)
    assert track.color == "#abcdef"
    assert track.region == TrackRegion("chr1", 100, 200)
    assert track.options == {"sample": "tumour"}
    with pytest.raises(TypeError):
        track.options["sample"] = "normal"


def test_track_canvas_exposes_stable_primitives_without_axes_access():
    fig, ax = plt.subplots()
    theme = load_config()
    canvas = TrackCanvas(
        ax, TrackRegion("chr1", 0, 100), "#123456",
        theme["visual_colors"],
    )

    canvas.set_y_limits(0, 1)
    canvas.set_y_ticks([0, 0.5, 1], ["0", ".5", "1"])
    canvas.horizontal_line(0.5)
    canvas.line([0, 100], [0.2, 0.8])
    canvas.step([0, 50, 100], [0.1, 0.9, 0.4])
    canvas.fill_between([0, 100], [0.2, 0.8])
    canvas.bars([20, 40], [0.3, 0.6], width=5)
    canvas.scatter([60], [0.7])
    canvas.rectangle(70, 0.1, 10, 0.3)
    canvas.span(82, 88)
    canvas.text(5, 0.9, "plugin")

    assert len(ax.lines) >= 3
    assert len(ax.patches) >= 2
    assert canvas.pixels_per_base > 0
    plt.close(fig)


def test_direct_module_plugin_loading_and_cli_spec_parsing():
    plugin = load_track_plugin("examples.gc_content_plugin:GCContentPlugin")
    assert plugin.name == "gc-content"

    source = build_plugin_track_sources([[
        "examples.gc_content_plugin:GCContentPlugin",
        "fasta=reference.fa", "window=20", "track_label=GC", "track_height=1.2",
        "track_color=#112233",
    ]])[0]

    assert source.label == "GC"
    assert source.height_in == pytest.approx(1.2)
    assert source.color == "#112233"
    assert source.options == {"fasta": "reference.fa", "window": "20"}


def test_installed_entry_point_discovery(monkeypatch):
    class FakeEntryPoint:
        name = "recording-ep"
        value = "example_package:RecordingPlugin"

        def load(self):
            return RecordingPlugin

    monkeypatch.setattr(
        plugin_api, "_entry_points_for_group", lambda: [FakeEntryPoint()]
    )

    assert plugin_api.available_track_plugins() == {
        "recording-ep": "example_package:RecordingPlugin"
    }
    assert load_track_plugin("recording-ep").name == "recording"


def test_unknown_or_incompatible_plugin_is_rejected():
    class FuturePlugin(RecordingPlugin):
        api_version = 99

    with pytest.raises(TrackPluginError, match="API version 99"):
        PluginTrackSource(FuturePlugin())
    with pytest.raises(TrackPluginError, match="KEY=VALUE"):
        build_plugin_track_sources([[
            "examples.gc_content_plugin:GCContentPlugin", "not-an-option"
        ]])


def test_renderer_draws_plugin_track_inside_standard_layout(tmp_path):
    track = PluginTrackSource(RecordingPlugin()).fetch("chr1", 100, 200)
    output = tmp_path / "plugin.png"
    renderer = AlignmentRenderer(
        show_alignments=False, show_coverage=False, show_ideogram=False,
        show_legend=False,
    )

    renderer.render(
        rows=[], chrom="chr1", window_start=100, window_end=200,
        reference=None, out_path=str(output), genomic_tracks=[track],
    )

    assert output.is_file()
    assert output.stat().st_size > 0


def test_bam_snapshot_accepts_public_plugin_tracks_argument(tmp_path):
    bam = os.path.join(os.path.dirname(__file__), "test.bam")
    snapshot = BamSnapshot(
        bam=bam, chrom="chr9", start=101867480, end=101867620,
        output_dir=str(tmp_path), output_name="plugin_snapshot",
        plugin_tracks=[PluginTrackSource(RecordingPlugin())],
        show_ideogram=False, show_coverage=False, show_legend=False,
        display_mode="squish", max_alignment_depth=10,
    )

    snapshot.snap()

    assert os.path.isfile(snapshot.output_path)


def test_plugin_fetch_and_render_failures_include_plugin_context(tmp_path):
    class BrokenPlugin(RecordingPlugin):
        name = "broken"

        def fetch(self, region, options):
            raise RuntimeError("remote service unavailable")

    with pytest.raises(TrackPluginError, match="broken.*fetching"):
        PluginTrackSource(BrokenPlugin()).fetch("chr1", 0, 10)

    class BrokenRenderPlugin(RecordingPlugin):
        name = "broken-render"

        def render(self, canvas, payload, region, options):
            raise RuntimeError("invalid payload")

    track = PluginTrackSource(BrokenRenderPlugin()).fetch("chr1", 0, 10)
    renderer = AlignmentRenderer(
        show_alignments=False, show_coverage=False, show_ideogram=False,
        show_legend=False,
    )
    with pytest.raises(TrackPluginError, match="broken-render.*rendering"):
        renderer.render(
            rows=[], chrom="chr1", window_start=0, window_end=10,
            reference=None, out_path=str(tmp_path / "broken.png"),
            genomic_tracks=[track],
        )
