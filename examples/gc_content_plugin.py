"""Minimal third-party-style GC-content track used by the demo gallery.

An external package can register this class under the
``locus_snap.track_plugins.v1`` entry-point group.  From a source checkout it
can also be loaded directly as ``examples.gc_content_plugin:GCContentPlugin``.
"""
from typing import Mapping

import pysam

from locus_snap import TrackCanvas, TrackPlugin, TrackRegion


class GCContentPlugin(TrackPlugin):
    name = "gc-content"
    default_label = "GC content plugin"
    default_height_in = 0.90
    default_color = "#188977"

    def fetch(self, region: TrackRegion, options: Mapping[str, str]):
        allowed = {"fasta", "window", "fill_alpha"}
        unknown = sorted(set(options) - allowed)
        if unknown:
            raise ValueError(f"unknown option(s): {', '.join(unknown)}")
        fasta = options.get("fasta")
        if not fasta:
            raise ValueError("the fasta=/path/reference.fa option is required")
        try:
            window = int(options.get("window", "25"))
        except ValueError as exc:
            raise ValueError("window must be an integer") from exc
        if window < 1:
            raise ValueError("window must be at least one base")

        with pysam.FastaFile(fasta) as reference:
            sequence = reference.fetch(region.chrom, region.start, region.end).upper()
        positions = []
        fractions = []
        for offset in range(0, len(sequence), window):
            chunk = sequence[offset:offset + window]
            called = sum(base in "ACGT" for base in chunk)
            gc = sum(base in "GC" for base in chunk)
            positions.append(region.start + offset + len(chunk) / 2)
            fractions.append(gc / called if called else 0.0)
        return {"positions": positions, "fractions": fractions}

    def render(
        self, canvas: TrackCanvas, payload, region: TrackRegion,
        options: Mapping[str, str],
    ) -> None:
        del region
        try:
            fill_alpha = float(options.get("fill_alpha", "0.24"))
        except ValueError as exc:
            raise ValueError("fill_alpha must be numeric") from exc
        if not 0 <= fill_alpha <= 1:
            raise ValueError("fill_alpha must be between zero and one")
        canvas.set_y_limits(0, 1)
        canvas.set_y_ticks([0, 0.5, 1], ["0", ".5", "1"])
        canvas.horizontal_line(0.5, linestyle="--", alpha=0.55)
        canvas.fill_between(
            payload["positions"], payload["fractions"],
            alpha=fill_alpha,
        )
        canvas.line(payload["positions"], payload["fractions"], linewidth=1.0)
