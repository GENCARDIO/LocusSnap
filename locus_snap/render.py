"""Matplotlib renderer.

Replaces the old approach of shelling out to `samtools tview`, capturing its
text table, and re-parsing that text with string splits. Here every read is
drawn from its own parsed CIGAR blocks, so insertions/deletions/soft-clips/
mismatches are geometrically exact instead of guessed from column spacing.
"""
from __future__ import annotations

from dataclasses import dataclass
from math import ceil, floor, log10, sqrt
from operator import attrgetter
from typing import Any, Dict, List, Optional, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import is_color_like, to_rgb
from matplotlib.lines import Line2D
from matplotlib.patches import ConnectionPatch, Patch, PathPatch, Polygon, Rectangle
from matplotlib.path import Path
from matplotlib.ticker import MaxNLocator

from locus_snap.annotations import (
    BAF_TRACK_FORMATS,
    CNV_TRACK_FORMATS,
    HIC_LOOP_TRACK_FORMATS,
    HIC_TRACK_FORMATS,
    PEAK_TRACK_FORMATS,
    SIGNAL_TRACK_FORMATS,
    TAD_TRACK_FORMATS,
    AnnotationItem,
    LoadedAnnotationTrack,
)
from locus_snap.config import (
    DEFAULT_CHROMOSOME_COLORS,
    DEFAULT_CHROMOSOME_PALETTE,
    DEFAULT_HAPLOTYPE_COLORS,
    DEFAULT_MODIFICATION_COLORS,
    DEFAULT_TAG_COLORS,
    DEFAULT_VISUAL_COLORS,
    load_config,
)
from locus_snap.cytobands import Cytoband
from locus_snap.read_model import AlignedRead, BaseModification
from locus_snap.reference import ReferenceWindow
from locus_snap.rna import (
    CANONICAL_SPLICE_MOTIFS,
    RNA_STRANDNESS_MODES,
    annotated_junctions,
    collect_fusion_evidence,
    collect_splice_junctions,
    gene_label_at,
    is_annotated_junction,
    splice_motif,
)
from locus_snap.track_plugin import (
    LoadedPluginTrack, TrackCanvas, TrackPluginError,
)

DEFAULT_COVERAGE_VAF_THRESHOLD = 0.20
DEFAULT_MAX_REFERENCE_SPAN = 250
MAX_EXPLICIT_MATE_CHROMOSOMES = 2
MAX_EXPLICIT_TAG_VALUES = 8
GRID_MODES = ("none", "major", "major_minor", "bands")
TITLE_ALIGNMENTS = ("left", "center", "right")

# Discordant-pair legend labels, IGV-equivalent: same categories/roles IGV's
# "color by insert size and pair orientation" mode uses (red = long insert,
# blue = short insert, a blue family for same-strand pairs, green for everted,
# per-chromosome hue for inter-chromosomal mates).
PAIR_CATEGORY_LABELS = {
    "large_insert": "Large insert",
    "small_insert": "Small insert",
    "ff": "FF (same strand)",
    "rr": "RR (same strand)",
    "everted": "Reverted (RF)",
    "interchrom": "Inter-chromosomal",
}


@dataclass
class SnvEvidence:
    count: int = 0
    forward: int = 0
    reverse: int = 0
    base_quality_sum: int = 0
    mapq_sum: int = 0


@dataclass
class ModificationEvidence:
    modified: int = 0
    canonical_depth: int = 0

    @property
    def fraction(self) -> float:
        return self.modified / self.canonical_depth if self.canonical_depth else 0.0


@dataclass(frozen=True)
class HighlightRegion:
    """A zero-based, half-open genomic interval highlighted across tracks."""
    chrom: str
    start: int
    end: int


def compute_feature_density(
    items: List[Any], window_start: int, window_end: int, bin_count: int
) -> Tuple[List[float], List[int]]:
    """Count interval overlaps in fixed-width bins across the visible window."""
    span = window_end - window_start
    if span <= 0 or bin_count <= 0:
        return [], []
    bin_count = min(bin_count, span)
    differences = [0] * (bin_count + 1)
    for item in items:
        lo = max(item.start, window_start)
        hi = min(item.end, window_end)
        if lo >= hi:
            continue
        first_bin = int((lo - window_start) * bin_count / span)
        final_bin = ceil((hi - window_start) * bin_count / span)
        first_bin = min(max(first_bin, 0), bin_count - 1)
        final_bin = min(max(final_bin, first_bin + 1), bin_count)
        differences[first_bin] += 1
        differences[final_bin] -= 1

    centers = []
    densities = []
    running = 0
    bin_width = span / bin_count
    for index in range(bin_count):
        running += differences[index]
        centers.append(window_start + (index + 0.5) * bin_width)
        densities.append(running)
    return centers, densities


def nice_scale_length(window_span: int) -> int:
    """Choose a UCSC-like 1/2/5 ruler length that fits most of the window."""
    if window_span <= 0:
        return 0
    target = max(window_span * 0.8, 1)
    magnitude = 10 ** floor(log10(target))
    candidates = [magnitude * multiplier for multiplier in (1, 2, 5, 10)]
    selected = candidates[0]
    for candidate in candidates:
        if candidate <= target:
            selected = candidate
    return min(max(int(selected), 1), window_span)


def format_scale_length(length: int) -> str:
    """Format an integer base-pair ruler length without scientific notation."""
    if length >= 1_000_000_000:
        return f"{length / 1_000_000_000:g} Gb"
    if length >= 1_000_000:
        return f"{length / 1_000_000:g} Mb"
    if length >= 1_000:
        return f"{length / 1_000:g} kb"
    return f"{length} bp"


def chrom_color(
    chrom: Optional[str], palette: Optional[List[str]] = None,
    *, colors: Optional[Dict[str, str]] = None,
) -> str:
    """Return IGV's mate-chromosome colour, with a stable contig fallback."""
    colors = DEFAULT_CHROMOSOME_COLORS if colors is None else colors
    palette = palette or DEFAULT_CHROMOSOME_PALETTE
    if not chrom:
        return DEFAULT_VISUAL_COLORS["axis"]
    label = str(chrom)
    normalized = label if label.startswith("chr") else f"chr{label}"
    if normalized in colors:
        return colors[normalized]
    h = 0
    for ch in normalized:
        h = (h * 31 + ord(ch)) & 0xFFFFFFFF
    return palette[h % len(palette)]


def categorical_color(
    category: Optional[str], colors: Dict[str, str],
    palette: Optional[List[str]] = None,
) -> str:
    palette = palette or DEFAULT_CHROMOSOME_PALETTE
    label = str(category) if category is not None else "untagged"
    if label in colors:
        return colors[label]
    if label.isdigit():
        return palette[(int(label) - 1) % len(palette)]
    value = 0
    for character in label:
        value = (value * 31 + ord(character)) & 0xFFFFFFFF
    return palette[value % len(palette)]


def haplotype_color(
    haplotype: Optional[str], colors: Optional[Dict[str, str]] = None,
    palette: Optional[List[str]] = None,
) -> str:
    return categorical_color(
        haplotype, colors or DEFAULT_HAPLOTYPE_COLORS, palette
    )


def tag_color(
    value: Optional[str], colors: Optional[Dict[str, str]] = None,
    palette: Optional[List[str]] = None,
) -> str:
    return categorical_color(value, colors or DEFAULT_TAG_COLORS, palette)


def modification_color(
    label: str, colors: Optional[Dict[str, str]] = None,
    palette: Optional[List[str]] = None,
) -> str:
    selected = colors or DEFAULT_MODIFICATION_COLORS
    if label in selected:
        return selected[label]
    if "other" in selected and len(selected) == 1:
        return selected["other"]
    return categorical_color(label, selected, palette)


def modification_matches(modification: BaseModification, codes: List[str]) -> bool:
    if not codes:
        return True
    accepted = {str(value).lower() for value in codes}
    aliases = {
        str(modification.code).lower(),
        modification.label.lower(),
        f"{modification.canonical_base}+{modification.code}".lower(),
    }
    return bool(accepted & aliases)


def has_base_modifications(
    reads: List[AlignedRead], start: int, end: int, codes: List[str],
) -> bool:
    return any(
        start <= modification.ref_position < end
        and modification_matches(modification, codes)
        for read in reads
        for modification in read.base_modifications
    )


def compute_modification_evidence(
    reads: List[AlignedRead], start: int, end: int,
    minimum_probability: float = 0.5, codes: Optional[List[str]] = None,
) -> Dict[str, Dict[int, ModificationEvidence]]:
    """Estimate modified/canonical depth from MM/ML calls.

    Calls at or above ``minimum_probability`` contribute to modified depth;
    every read carrying the canonical base at that genomic position contributes
    to the denominator. The result is intentionally sparse and contains only
    positions represented by at least one MM/ML call.
    """
    selected_codes = list(codes or [])
    evidence: Dict[str, Dict[int, ModificationEvidence]] = {}
    canonical_by_site: Dict[Tuple[str, int], set] = {}
    seen_calls = set()
    for read_index, read in enumerate(reads):
        for modification in read.base_modifications:
            if not start <= modification.ref_position < end:
                continue
            if not modification_matches(modification, selected_codes):
                continue
            label = modification.label
            position = modification.ref_position
            canonical_by_site.setdefault((label, position), set()).add(
                modification.aligned_base
            )
            item = evidence.setdefault(label, {}).setdefault(
                position, ModificationEvidence()
            )
            call_key = (read_index, label, position)
            if call_key in seen_calls:
                continue
            seen_calls.add(call_key)
            if (
                modification.probability is not None
                and modification.probability >= minimum_probability
            ):
                item.modified += 1

    for (label, position), canonical_bases in canonical_by_site.items():
        item = evidence[label][position]
        item.canonical_depth = sum(
            1 for read in reads if read.base_at(position) in canonical_bases
        )
    return evidence


def compute_coverage(reads: List[AlignedRead], start: int, end: int) -> List[int]:
    """Per-base depth across [start, end), counting only reference-consuming
    match bases (a deleted base is not "covered", matching `samtools depth`)."""
    span = max(0, end - start)
    changes = [0] * (span + 1)
    for read in reads:
        for blk in read.blocks:
            if blk.op not in ("M", "=", "X"):
                continue
            lo = max(blk.ref_pos, start)
            hi = min(blk.ref_pos + blk.length, end)
            if lo < hi:
                changes[lo - start] += 1
                changes[hi - start] -= 1
    depth = []
    running = 0
    for change in changes[:-1]:
        running += change
        depth.append(running)
    return depth


def compute_binned_coverage(
    reads: List[AlignedRead], start: int, end: int, bin_count: int
) -> Tuple[List[float], List[float]]:
    """Return bin edges and mean depth without allocating one value per base."""
    span = max(0, end - start)
    if span == 0 or bin_count <= 0:
        return [], []
    bin_count = min(bin_count, span)
    bin_width = span / bin_count
    depth_sums = [0.0] * bin_count
    for read in reads:
        for block in read.blocks:
            if block.op not in ("M", "=", "X"):
                continue
            lo = max(block.ref_pos, start)
            hi = min(block.ref_pos + block.length, end)
            if lo >= hi:
                continue
            first_bin = min(int((lo - start) / bin_width), bin_count - 1)
            final_bin = min(ceil((hi - start) / bin_width), bin_count)
            for bin_index in range(first_bin, final_bin):
                bin_start = start + bin_index * bin_width
                bin_end = bin_start + bin_width
                overlap = max(0.0, min(hi, bin_end) - max(lo, bin_start))
                depth_sums[bin_index] += overlap

    edges = []
    depths = []
    for bin_index in range(bin_count + 1):
        edges.append(start + bin_index * bin_width)
    for depth_sum in depth_sums:
        depths.append(depth_sum / bin_width)
    return edges, depths


def compute_snv_counts(
    reads: List[AlignedRead], start: int, end: int
) -> Dict[int, Dict[str, int]]:
    """Count reference-backed A/C/G/T mismatches at each covered position."""
    counts: Dict[int, Dict[str, int]] = {}
    for read in reads:
        for position, base in read.mismatches:
            if start <= position < end and base in "ACGT":
                position_counts = counts.setdefault(position, {})
                position_counts[base] = position_counts.get(base, 0) + 1
    return counts


def compute_snv_evidence(
    reads: List[AlignedRead], start: int, end: int,
    min_baseq: int = 0, min_mapq: int = 0,
) -> Tuple[List[int], Dict[int, Dict[str, SnvEvidence]]]:
    """Return quality-filtered nucleotide depth and alternate-SNV evidence."""
    depth = [0] * max(0, end - start)
    evidence: Dict[int, Dict[str, SnvEvidence]] = {}
    for read in reads:
        if getattr(read, "mapq", 0) < min_mapq:
            continue
        sequence = getattr(read, "query_sequence", "") or ""
        qualities = getattr(read, "query_qualities", []) or []
        for block in read.blocks:
            if block.op not in ("M", "=", "X"):
                continue
            lo = max(block.ref_pos, start)
            hi = min(block.ref_pos + block.length, end)
            for position in range(lo, hi):
                query_index = block.query_pos + position - block.ref_pos
                base = sequence[query_index].upper() if query_index < len(sequence) else "N"
                base_quality = qualities[query_index] if query_index < len(qualities) else 0
                if base in "ACGT" and base_quality >= min_baseq:
                    depth[position - start] += 1

        details = getattr(read, "mismatch_details", None)
        if details is None:
            details = []
            for position, base in getattr(read, "mismatches", []):
                details.append((position, base, 0))
        for position, base, base_quality in details:
            if not start <= position < end or base not in "ACGT" or base_quality < min_baseq:
                continue
            base_evidence = evidence.setdefault(position, {}).setdefault(base, SnvEvidence())
            base_evidence.count += 1
            if getattr(read, "is_reverse", False):
                base_evidence.reverse += 1
            else:
                base_evidence.forward += 1
            base_evidence.base_quality_sum += base_quality
            base_evidence.mapq_sum += getattr(read, "mapq", 0)
    return depth, evidence


def compute_sparse_snv_evidence(
    reads: List[AlignedRead], start: int, end: int,
    min_baseq: int = 0, min_mapq: int = 0,
) -> Tuple[Dict[int, int], Dict[int, Dict[str, SnvEvidence]]]:
    """Compute SNV evidence without a window-sized depth array.

    Only positions carrying qualifying mismatch evidence need a depth
    denominator, which keeps wide, mostly empty windows memory-bounded.
    """
    evidence: Dict[int, Dict[str, SnvEvidence]] = {}
    for read in reads:
        if getattr(read, "mapq", 0) < min_mapq:
            continue
        details = getattr(read, "mismatch_details", None) or []
        for position, base, base_quality in details:
            if not start <= position < end or base not in "ACGT" or base_quality < min_baseq:
                continue
            base_evidence = evidence.setdefault(position, {}).setdefault(base, SnvEvidence())
            base_evidence.count += 1
            if getattr(read, "is_reverse", False):
                base_evidence.reverse += 1
            else:
                base_evidence.forward += 1
            base_evidence.base_quality_sum += base_quality
            base_evidence.mapq_sum += getattr(read, "mapq", 0)

    if not evidence:
        return {}, {}
    candidate_positions = set(evidence)
    depth = {}
    for position in candidate_positions:
        depth[position] = 0
    for read in reads:
        if getattr(read, "mapq", 0) < min_mapq:
            continue
        sequence = getattr(read, "query_sequence", "") or ""
        qualities = getattr(read, "query_qualities", []) or []
        for block in read.blocks:
            if block.op not in ("M", "=", "X"):
                continue
            lo = max(block.ref_pos, start)
            hi = min(block.ref_pos + block.length, end)
            for position in range(lo, hi):
                if position not in candidate_positions:
                    continue
                query_index = block.query_pos + position - block.ref_pos
                base = sequence[query_index].upper() if query_index < len(sequence) else "N"
                base_quality = qualities[query_index] if query_index < len(qualities) else 0
                if base in "ACGT" and base_quality >= min_baseq:
                    depth[position] += 1
    return depth, evidence


def compute_splice_junctions(
    reads: List[AlignedRead], strand_mode: str = "combined"
) -> Dict[Tuple[int, int, str], int]:
    """Backward-compatible count view over richer junction evidence."""
    return {
        key: evidence.count
        for key, evidence in collect_splice_junctions(
            reads, strand_mode=strand_mode
        ).items()
    }


def nice_tick_positions(start: int, end: int, target: int = 8) -> List[int]:
    locator = MaxNLocator(nbins=target, steps=[1, 2, 5, 10])
    ticks = []
    for value in locator.tick_values(start, end):
        if start <= value <= end:
            ticks.append(int(value))
    return ticks


def genomic_tick_labels(
    ticks: List[int], window_length: int,
) -> Tuple[List[str], str]:
    """Format positions using the unit appropriate for the window span."""
    span = max(window_length, 1)
    if span >= 1_000_000_000:
        scale, unit = 1_000_000_000, "Gb"
    elif span >= 1_000_000:
        scale, unit = 1_000_000, "Mb"
    elif span >= 1_000:
        scale, unit = 1_000, "kb"
    else:
        scale, unit = 1, "bp"

    minimum_step = None
    for left, right in zip(ticks, ticks[1:]):
        step = abs(right - left)
        if step and (minimum_step is None or step < minimum_step):
            minimum_step = step
    scaled_step = (minimum_step or scale) / scale
    decimals = 0
    while scaled_step < 1 and decimals < 9:
        scaled_step *= 10
        decimals += 1

    labels = []
    for tick in ticks:
        labels.append(f"{tick / scale:,.{decimals}f}")
    while len(set(labels)) < len(labels) and decimals < 9:
        decimals += 1
        labels = []
        for tick in ticks:
            labels.append(f"{tick / scale:,.{decimals}f}")
    return labels, unit


def apply_genomic_axis(
    ax, ticks: List[int], window_start: int, window_end: int,
    label_size: float = 9, color: str = "#0b0b0b",
) -> None:
    """Apply fixed genomic labels and suppress Matplotlib's scientific offset."""
    labels, unit = genomic_tick_labels(ticks, window_end - window_start + 1)
    ax.set_xticks(ticks)
    ax.set_xticklabels(labels)
    ax.set_xlabel(
        f"Position ({unit})", fontsize=max(label_size - 1, 6),
        color=color, labelpad=4,
    )
    ax.get_xaxis().get_offset_text().set_visible(False)


def left_margin_fraction(
    fig_width: float, genomic_tracks: List[Any]
) -> float:
    """Reserve enough physical space for annotation labels outside the axes."""
    longest_label = 0
    for track in genomic_tracks:
        if isinstance(track, LoadedPluginTrack):
            # Plugin titles sit inside their axes so their optional numeric
            # y ticks retain the outside margin without collisions.
            continue
        if len(track.label) > longest_label:
            longest_label = len(track.label)
    margin_in = max(0.70, 0.25 + longest_label * 0.065)
    return min(margin_in / fig_width, 0.22)


def annotation_row_count(track: Any) -> int:
    """Return a layout row count for built-in or plugin-backed tracks."""
    if isinstance(track, LoadedPluginTrack):
        return 1
    return max(len(track.rows), 1)


def ellipsize(text: str, max_chars: int) -> str:
    """Fit a label to an approximate character budget without spilling."""
    if max_chars <= 0:
        return ""
    if len(text) <= max_chars:
        return text
    if max_chars == 1:
        return "…"
    return text[:max_chars - 1].rstrip() + "…"


class AlignmentRenderer:
    def __init__(
        self,
        fig_width: float = 14.0,
        dpi: int = 150,
        show_alignments: bool = True,
        show_legend: bool = True,
        show_coverage: bool = True,
        annotate_gap: bool = True,
        max_mismatch_render_span: int = 5000,
        pair_colors: bool = True,
        shade_by_mapq: bool = True,
        mapq_cap: int = 60,
        alignment_colors: Optional[Dict[str, Optional[str]]] = None,
        display_mode: str = "expand",
        show_ideogram: bool = True,
        max_reference_span: int = DEFAULT_MAX_REFERENCE_SPAN,
        view_as_pairs: bool = False,
        coverage_vaf_threshold: float = DEFAULT_COVERAGE_VAF_THRESHOLD,
        min_baseq: int = 0,
        min_variant_mapq: int = 0,
        show_variant_counts: bool = False,
        show_indel_lengths: bool = False,
        haplotype_view: str = "none",
        read_tag: Optional[str] = None,
        tag_view: str = "none",
        tag_label: Optional[str] = None,
        tag_colors: Optional[Dict[str, str]] = None,
        visual_config: Optional[Dict[str, Any]] = None,
        sort_base_position: Optional[int] = None,
        sort_reference_base: Optional[str] = None,
        show_center_guide: bool = False,
        show_sashimi: bool = False,
        min_junction_reads: int = 1,
        sashimi_strand: str = "combined",
        min_junction_anchor: int = 0,
        rna_strandness: str = "alignment",
        junction_labels: str = "count",
        show_fusions: bool = False,
        min_fusion_reads: int = 2,
        fusion_breakpoint_tolerance: int = 10,
        fusion_min_distance: int = 100_000,
        min_fusion_mapq: int = 20,
        grid_mode: str = "major",
        highlight_regions: Optional[List[HighlightRegion]] = None,
        highlight_color: str = "#ffd54f",
        highlight_alpha: float = 0.20,
        title_align: str = "left",
        long_read_mode: bool = False,
        show_base_modifications: bool = False,
        modification_codes: Optional[List[str]] = None,
        min_mod_probability: float = 0.5,
        molecule_mode: bool = False,
    ):
        theme = visual_config or load_config()
        self.base_colors = dict(theme["base_colors"])
        self.track_colors = dict(theme["track_colors"])
        self.visual_colors = dict(theme["visual_colors"])
        self.haplotype_colors = dict(theme["haplotype_colors"])
        self.tag_colors = dict(theme["tag_colors"])
        self.tag_colors.update(tag_colors or {})
        for category, color in self.tag_colors.items():
            if not is_color_like(color):
                raise ValueError(f"Invalid color for tag value {category!r}: {color!r}.")
        self.cytoband_colors = dict(theme["cytoband_colors"])
        self.chromosome_colors = dict(theme["chromosome_colors"])
        self.chromosome_palette = list(theme["chromosome_palette"])
        self.long_read_colors = dict(theme["long_read_colors"])
        self.modification_colors = dict(theme["modification_colors"])
        self.molecule_colors = dict(theme["molecule_colors"])
        self.styles = dict(theme["styles"])
        self.sort_base_position = sort_base_position
        self.sort_reference_base = sort_reference_base
        self.active_sort_base_position = sort_base_position
        self.active_sort_reference_base = sort_reference_base
        self.show_center_guide = show_center_guide
        self.show_sashimi = show_sashimi
        self.show_fusions = show_fusions
        self.show_rna_evidence = show_sashimi or show_fusions
        if min_junction_reads < 1:
            raise ValueError("Minimum junction-read support must be at least one.")
        self.min_junction_reads = min_junction_reads
        if sashimi_strand not in ("combined", "split"):
            raise ValueError("Sashimi strand mode must be combined or split.")
        self.sashimi_strand = sashimi_strand
        if min_junction_anchor < 0:
            raise ValueError("Minimum junction anchor cannot be negative.")
        self.min_junction_anchor = min_junction_anchor
        if rna_strandness not in RNA_STRANDNESS_MODES:
            raise ValueError(
                f"RNA strandness must be one of: {', '.join(RNA_STRANDNESS_MODES)}."
            )
        self.rna_strandness = rna_strandness
        if junction_labels not in ("count", "status", "full"):
            raise ValueError("Junction labels must be count, status, or full.")
        self.junction_labels = junction_labels
        if min_fusion_reads < 1:
            raise ValueError("Minimum fusion-read support must be at least one.")
        if fusion_breakpoint_tolerance < 0 or fusion_min_distance < 0:
            raise ValueError("Fusion breakpoint tolerance and distance cannot be negative.")
        if min_fusion_mapq < 0:
            raise ValueError("Minimum fusion MAPQ cannot be negative.")
        self.min_fusion_reads = min_fusion_reads
        self.fusion_breakpoint_tolerance = fusion_breakpoint_tolerance
        self.fusion_min_distance = fusion_min_distance
        self.min_fusion_mapq = min_fusion_mapq
        self.fig_width = fig_width
        self.dpi = dpi
        self.show_alignments = show_alignments
        self.show_legend = show_legend
        if grid_mode not in GRID_MODES:
            raise ValueError(
                f"Unknown grid mode '{grid_mode}'. Choose {', '.join(GRID_MODES)}."
            )
        self.grid_mode = grid_mode
        if not is_color_like(highlight_color):
            raise ValueError(f"Invalid highlight color: {highlight_color!r}.")
        if (
            isinstance(highlight_alpha, bool)
            or not isinstance(highlight_alpha, (int, float))
            or not 0 < highlight_alpha <= 1
        ):
            raise ValueError("Highlight alpha must be greater than 0 and at most 1.")
        self.highlight_regions = list(highlight_regions or [])
        for region in self.highlight_regions:
            if region.end <= region.start:
                raise ValueError("Highlight region end must be greater than its start.")
        self.highlight_color = highlight_color
        self.highlight_alpha = float(highlight_alpha)
        if title_align not in TITLE_ALIGNMENTS:
            raise ValueError(
                f"Unknown title alignment '{title_align}'. "
                f"Choose {', '.join(TITLE_ALIGNMENTS)}."
            )
        self.title_align = title_align
        self.long_read_mode = long_read_mode
        self.show_base_modifications = show_base_modifications
        self.modification_codes = list(modification_codes or [])
        if not 0 <= min_mod_probability <= 1:
            raise ValueError("Minimum modification probability must be between 0 and 1.")
        self.min_mod_probability = float(min_mod_probability)
        self.molecule_mode = molecule_mode
        self.modification_labels_seen: set = set()
        self.show_coverage = show_coverage
        self.annotate_gap = annotate_gap
        self.max_mismatch_render_span = max_mismatch_render_span
        self.pair_colors = pair_colors and not long_read_mode and not molecule_mode
        self.shade_by_mapq = shade_by_mapq
        self.mapq_cap = mapq_cap
        self.alignment_colors = dict(theme["alignment_colors"])
        if alignment_colors:
            self.alignment_colors.update(alignment_colors)
        self.interchrom_mate_colors: Dict[str, str] = {}
        if display_mode not in ("collapse", "expand", "squish"):
            raise ValueError(
                f"Unknown display mode '{display_mode}'. Choose collapse, expand, or squish."
            )
        self.display_mode = display_mode
        self.show_ideogram = show_ideogram
        self.max_reference_span = max_reference_span
        self.view_as_pairs = view_as_pairs
        if not 0 <= coverage_vaf_threshold <= 1:
            raise ValueError("Coverage VAF threshold must be between 0 and 1.")
        self.coverage_vaf_threshold = coverage_vaf_threshold
        if min_baseq < 0 or min_variant_mapq < 0:
            raise ValueError("Variant base-quality and MAPQ filters cannot be negative.")
        self.min_baseq = min_baseq
        self.min_variant_mapq = min_variant_mapq
        self.show_variant_counts = show_variant_counts
        self.show_indel_lengths = show_indel_lengths
        if haplotype_view not in ("none", "color", "split"):
            raise ValueError("Haplotype view must be none, color, or split.")
        if tag_view not in ("none", "color", "split"):
            raise ValueError("Tag view must be none, color, or split.")
        if haplotype_view != "none" and tag_view != "none":
            raise ValueError("Haplotype and generic tag views cannot be active together.")
        if tag_view != "none" and not read_tag:
            raise ValueError("Tag colour/group view requires a BAM tag name.")
        self.haplotype_view = haplotype_view
        self.read_tag = read_tag
        self.tag_view = tag_view
        self.long_read_coloring = (
            long_read_mode and haplotype_view == "none" and tag_view == "none"
        )
        self.tag_label = tag_label or (f"Tag {read_tag}" if read_tag else "Tag")
        self.tag_value_colors: Dict[str, str] = {}
        self.tag_value_counts: Dict[str, int] = {}
        self.has_split_lanes = haplotype_view == "split" or tag_view == "split"
        self.row_height_in = (
            self.styles["squish_row_height_in"] if display_mode == "squish"
            else self.styles["row_height_in"]
        )
        self.row_margin = (
            self.styles["squish_row_margin"] if display_mode == "squish"
            else self.styles["row_margin"]
        )
        if self.long_read_coloring and show_base_modifications:
            legend_group_count = 5
        elif self.long_read_coloring or show_base_modifications:
            legend_group_count = 4
        elif molecule_mode or haplotype_view != "none" or tag_view != "none":
            legend_group_count = 4
        elif self.pair_colors:
            legend_group_count = 4
        else:
            legend_group_count = 3
        if fig_width >= 9:
            self.legend_height_in = 0.86
        elif fig_width >= 6:
            legend_rows = (legend_group_count + 1) // 2
            self.legend_height_in = legend_rows * 0.72
        else:
            self.legend_height_in = legend_group_count * 0.62
        legend_enabled = show_alignments and show_legend
        self.legend_bottom_in = 0.04 if legend_enabled else 0.0
        self.legend_plot_gap_in = 0.12 if legend_enabled else 0.0
        self.legend_tick_clearance_in = 0.30 if legend_enabled else 0.44
        if not legend_enabled:
            self.legend_height_in = 0.0
        self.legend_margin_in = (
            self.legend_bottom_in + self.legend_height_in + self.legend_plot_gap_in
            + self.legend_tick_clearance_in
        )

    def draw_background_grid(
        self,
        ax,
        major_ticks: List[int],
        window_start: int,
        window_end: int,
    ) -> None:
        """Draw a coordinate-aware background shared by genomic tracks.

        Minor lines and alternating bands are derived from the labelled major
        coordinates, so every track and panel uses the same genomic cadence.
        """
        if self.grid_mode == "none" or window_end <= window_start:
            return

        ticks = sorted({
            float(tick) for tick in major_ticks
            if window_start <= tick <= window_end
        })
        boundaries = sorted({float(window_start), *ticks, float(window_end)})

        if self.grid_mode == "bands":
            for index, (left, right) in enumerate(zip(boundaries, boundaries[1:])):
                if index % 2 == 0 and right > left:
                    ax.axvspan(
                        left, right,
                        facecolor=self.visual_colors["grid_band"],
                        edgecolor="none",
                        alpha=self.styles["grid_band_alpha"],
                        zorder=-1,
                    )

        if self.grid_mode == "major_minor" and len(boundaries) >= 2:
            divisions = max(1, int(round(self.styles["grid_minor_divisions"])))
            if divisions > 1:
                for left, right in zip(boundaries, boundaries[1:]):
                    step = (right - left) / divisions
                    for subdivision in range(1, divisions):
                        position = left + subdivision * step
                        ax.axvline(
                            position,
                            color=self.visual_colors["gridline_minor"],
                            alpha=self.styles["minor_grid_line_alpha"],
                            lw=self.styles["minor_grid_line_width"],
                            linestyle=self.styles["minor_grid_line_style"],
                            zorder=0,
                        )

        for tick in ticks:
            ax.axvline(
                tick,
                color=self.visual_colors["gridline"],
                alpha=self.styles["grid_line_alpha"],
                lw=self.styles["grid_line_width"],
                linestyle=self.styles["grid_line_style"],
                zorder=0,
            )

    def draw_highlights(
        self,
        ax,
        chrom: str,
        window_start: int,
        window_end: int,
    ) -> None:
        """Shade configured intervals behind the data in one genomic track."""
        normalized_chrom = chrom.lower()
        if normalized_chrom.startswith("chr"):
            normalized_chrom = normalized_chrom[3:]
        for region in self.highlight_regions:
            region_chrom = region.chrom.lower()
            if region_chrom.startswith("chr"):
                region_chrom = region_chrom[3:]
            if region_chrom != normalized_chrom:
                continue
            left = max(region.start, window_start)
            right = min(region.end, window_end)
            if right <= left:
                continue
            ax.axvspan(
                left, right,
                facecolor=self.highlight_color,
                edgecolor="none",
                alpha=self.highlight_alpha,
                zorder=0.25,
            )

    def figure_title_position(self) -> Tuple[float, str]:
        """Return the figure x-coordinate and anchor for the title block."""
        if self.title_align == "center":
            return 0.50, "center"
        if self.title_align == "right":
            return 0.99, "right"
        return 0.01, "left"

    def annotation_track_height(
        self, annotation: Any, row_count: Optional[int] = None,
    ) -> float:
        if isinstance(annotation, LoadedPluginTrack):
            return annotation.height_in
        """Resolve a custom override or the YAML default for an annotation track."""
        if annotation.height_in is not None:
            return annotation.height_in
        if (
            annotation.display_mode == "density"
            or annotation.kind in PEAK_TRACK_FORMATS
            or annotation.kind in SIGNAL_TRACK_FORMATS
        ):
            return self.styles["peak_track_height_in"]
        if annotation.kind in CNV_TRACK_FORMATS:
            return self.styles["cnv_track_height_in"]
        if annotation.kind in BAF_TRACK_FORMATS:
            return self.styles["baf_track_height_in"]
        if annotation.kind in HIC_TRACK_FORMATS:
            return self.styles["hic_track_height_in"]
        rows = row_count if row_count is not None else max(len(annotation.rows), 1)
        return max(rows, 1) * self.styles["annotation_row_height_in"]

    def read_style(self, read: AlignedRead):
        """(fill_color, alpha) for a read's main body: hue encodes pair
        discordance category (when enabled), alpha encodes mapping quality -
        low-MAPQ reads get a lighter/more washed-out fill, same idea as IGV's
        "shade by mapping quality"."""
        if self.haplotype_view in ("color", "split"):
            color = haplotype_color(
                getattr(read, "haplotype", None), self.haplotype_colors,
                self.chromosome_palette,
            )
        elif self.tag_view in ("color", "split"):
            raw_value = getattr(read, "tag_value", None)
            label = str(raw_value) if raw_value is not None else "untagged"
            color = tag_color(
                raw_value, self.tag_colors, self.chromosome_palette,
            )
            self.tag_value_colors[label] = color
            self.tag_value_counts[label] = self.tag_value_counts.get(label, 0) + 1
        elif self.long_read_coloring:
            if read.is_supplementary:
                color = self.long_read_colors["supplementary"]
            elif read.is_reverse:
                color = self.long_read_colors["reverse"]
            else:
                color = self.long_read_colors["forward"]
        elif self.molecule_mode:
            if getattr(read, "molecule_is_duplex", False):
                color = self.molecule_colors["duplex"]
            elif getattr(read, "molecule_family_size", 1) > 1:
                color = self.molecule_colors["consensus"]
            else:
                color = self.molecule_colors["singleton"]
        elif self.pair_colors and read.pair_category == "interchrom":
            color = self.alignment_colors["interchrom"] or chrom_color(
                read.mate_chrom, self.chromosome_palette,
                colors=self.chromosome_colors,
            )
            if read.mate_chrom:
                self.interchrom_mate_colors[str(read.mate_chrom)] = color
        elif self.pair_colors and read.pair_category == "large_insert":
            color = self.alignment_colors["large_insert"]
        elif self.pair_colors and read.pair_category == "small_insert":
            color = self.alignment_colors["small_insert"]
        elif self.pair_colors and read.pair_category == "ff":
            color = self.alignment_colors["ff"]
        elif self.pair_colors and read.pair_category == "rr":
            color = self.alignment_colors["rr"]
        elif self.pair_colors and read.pair_category == "everted":
            color = self.alignment_colors["everted"]
        else:
            color = self.alignment_colors["normal"]

        alpha = (
            self.styles["secondary_alignment_alpha"]
            if (read.is_secondary or read.is_duplicate)
            else self.styles["alignment_alpha"]
        )
        if self.shade_by_mapq and self.mapq_cap > 0:
            mapq_frac = min(max(read.mapq, 0), self.mapq_cap) / self.mapq_cap
            floor = self.styles["mapq_alpha_floor"]
            alpha *= floor + (1 - floor) * mapq_frac
        return color, alpha

    def render(
        self,
        rows: List[List[AlignedRead]],
        chrom: str,
        window_start: int,
        window_end: int,
        reference: Optional[ReferenceWindow],
        out_path: str,
        title: str = "",
        layout: str = "pack",
        dropped_reads: int = 0,
        downsampled_reads: int = 0,
        all_reads_for_coverage: Optional[List[AlignedRead]] = None,
        genomic_tracks: Optional[List[LoadedAnnotationTrack]] = None,
        contig_length: Optional[int] = None,
        cytobands: Optional[List[Cytoband]] = None,
        assembly_label: Optional[str] = None,
    ) -> None:
        span = window_end - window_start
        n_rows = max(len(rows), 1)
        show_ref_track = bool(
            reference and reference.available and
            self.max_reference_span > 0 and span <= self.max_reference_span
        )
        render_base_detail = span <= self.max_mismatch_render_span
        data_reads = all_reads_for_coverage
        if data_reads is None:
            data_reads = [read for row in rows for read in row]
        show_modification_track = bool(
            self.show_base_modifications
            and has_base_modifications(
                data_reads, window_start, window_end, self.modification_codes
            )
        )

        tracks = []
        ratios = []
        if self.show_ideogram and contig_length:
            tracks.append("ideogram")
            ratios.append(self.styles["ideogram_height_in"])
        if show_ref_track:
            tracks.append("reference")
            ratios.append(self.styles["reference_height_in"])
        genomic_tracks = genomic_tracks or []
        for index, annotation in enumerate(genomic_tracks):
            tracks.append(f"annotation_{index}")
            ratios.append(self.annotation_track_height(annotation))
        if self.show_coverage:
            tracks.append("coverage")
            ratios.append(self.styles["coverage_track_height_in"])
        if show_modification_track:
            tracks.append("modifications")
            ratios.append(self.styles["modification_track_height_in"])
        if self.show_rna_evidence:
            tracks.append("sashimi")
            ratios.append(self.styles["sashimi_track_height_in"])
        if self.show_alignments:
            tracks.append("alignments")
            ratios.append(max(n_rows * self.row_height_in, self.row_height_in))
        if not tracks:
            raise ValueError(
                "Nothing to render: enable alignments, coverage, the ideogram, or a genomic track."
            )

        top_margin_in = 0.72  # dedicated region-title and subtitle lanes
        bottom_margin_in = self.legend_margin_in
        fig_height = sum(ratios) + top_margin_in + bottom_margin_in
        fig, axes = plt.subplots(
            nrows=len(tracks),
            ncols=1,
            figsize=(self.fig_width, fig_height),
            dpi=self.dpi,
            gridspec_kw={"height_ratios": ratios, "hspace": 0.15},
            sharex=True,
        )
        if len(tracks) == 1:
            axes = [axes]
        ax_by_track = dict(zip(tracks, axes))

        for ax in axes:
            ax.set_xlim(window_start, window_end)
            for spine in ("top", "right", "left"):
                ax.spines[spine].set_visible(False)
            ax.spines["bottom"].set_visible(False)
            ax.tick_params(
                left=False, labelleft=False, bottom=False, top=False,
                labelbottom=False, labeltop=False,
            )

        tick_positions = nice_tick_positions(window_start, window_end)

        title_x, title_ha = self.figure_title_position()
        fig.text(
            title_x, 1 - 0.06 / fig_height,
            f"{chrom}:{window_start + 1:,}-{window_end:,} ({span:,} bp)",
            fontsize=10.5, color=self.visual_colors["primary_text"], fontweight="bold", va="top", ha=title_ha,
        )
        subtitle = title
        if dropped_reads:
            subtitle = (subtitle + " -- " if subtitle else "") + (
                f"{dropped_reads} lower-priority read(s) not shown (--max_rows)"
            )
        if downsampled_reads:
            subtitle = (subtitle + " -- " if subtitle else "") + (
                f"{downsampled_reads} alignment(s) downsampled"
            )
        if subtitle:
            fig.text(
                title_x, 1 - 0.34 / fig_height,
                ellipsize(subtitle, max(30, int(self.fig_width * 15))),
                fontsize=8.5, color=self.visual_colors["secondary_text"], va="top", ha=title_ha,
            )

        if "ideogram" in ax_by_track:
            self.draw_ideogram(
                ax_by_track["ideogram"], chrom, window_start, window_end, contig_length,
                cytobands,
            )

        # --- reference -----------------------------------------------------
        if show_ref_track:
            self.draw_reference_track(
                ax_by_track["reference"], reference, window_start, window_end,
                available_width_in=self.fig_width,
            )

        for index, annotation in enumerate(genomic_tracks):
            self.draw_annotation_track(
                ax_by_track[f"annotation_{index}"], annotation, window_start, window_end
            )

        # --- coverage --------------------------------------------------
        if self.show_coverage:
            cov_ax = ax_by_track["coverage"]
            self.draw_coverage_track(cov_ax, data_reads, window_start, window_end)

        if show_modification_track:
            self.draw_modification_track(
                ax_by_track["modifications"], data_reads, window_start, window_end
            )

        if self.show_rna_evidence:
            sashimi_reads = all_reads_for_coverage
            if sashimi_reads is None:
                sashimi_reads = []
                for row in rows:
                    sashimi_reads.extend(row)
            self.draw_sashimi_track(
                ax_by_track["sashimi"], sashimi_reads, window_start, window_end,
                chrom=chrom, reference=reference, genomic_tracks=genomic_tracks,
            )

        # --- alignments and genomic axis --------------------------------
        axis_ax = ax_by_track[tracks[-1]]
        if self.show_alignments:
            aln_ax = ax_by_track["alignments"]
            aln_ax.set_ylim(n_rows, 0)
            self.draw_haplotype_lanes(aln_ax, rows)
            for row_idx, row in enumerate(rows):
                y0 = row_idx + self.row_margin
                h = 1 - 2 * self.row_margin
                self.draw_alignment_row(
                    aln_ax, row, y0, h, render_base_detail, layout
                )

            if not rows:
                aln_ax.text(
                    0.5, 0.5, "No alignments in this region", transform=aln_ax.transAxes,
                    ha="center", va="center", fontsize=10, color=self.visual_colors["secondary_text"],
                )
        apply_genomic_axis(
            axis_ax, tick_positions, window_start, window_end, label_size=9,
            color=self.visual_colors["primary_text"],
        )
        axis_ax.tick_params(
            bottom=True, labelbottom=True, labelsize=9,
            length=3, colors=self.visual_colors["primary_text"],
        )
        for track in tracks:
            if track != "ideogram":
                self.draw_background_grid(
                    ax_by_track[track], tick_positions, window_start, window_end
                )
                self.draw_highlights(
                    ax_by_track[track], chrom, window_start, window_end
                )

        for track, ax in ax_by_track.items():
            if track != "ideogram":
                self.draw_center_guide(ax, window_start, window_end)

        # --- legend -----------------------------------
        plot_left = left_margin_fraction(self.fig_width, genomic_tracks)
        if self.has_split_lanes:
            plot_left = max(plot_left, min(1.15 / self.fig_width, 0.25))
        plot_right = 0.92
        fig.subplots_adjust(left=plot_left, right=plot_right,
                            top=1 - top_margin_in / fig_height,
                            bottom=bottom_margin_in / fig_height)
        self.draw_scale_bar(
            fig, fig_height, plot_left, plot_right, span, assembly_label,
            offset_from_top_in=0.59,
        )
        if self.show_alignments and self.show_legend:
            legends = self.draw_legends(fig, fig_height, plot_left, plot_right)
            self.separate_legend_from_plots(fig, axes, legends)
        fig.savefig(out_path, dpi=self.dpi)
        plt.close(fig)

    def draw_scale_bar(
        self,
        fig,
        fig_height: float,
        plot_left: float,
        plot_right: float,
        window_span: int,
        assembly_label: Optional[str] = None,
        offset_from_top_in: float = 0.59,
    ) -> None:
        """Draw a UCSC-style labelled genomic ruler in the figure margin."""
        scale_length = nice_scale_length(window_span)
        if scale_length <= 0 or plot_right <= plot_left:
            return
        line_start = plot_left
        line_end = plot_left + (plot_right - plot_left) * scale_length / window_span
        line_y = 1 - offset_from_top_in / fig_height
        cap_half_height = 0.055 / fig_height
        line_color = self.visual_colors["primary_text"]
        line_width = 0.8
        fig.lines.extend([
            Line2D(
                [line_start, line_end], [line_y, line_y],
                transform=fig.transFigure, color=line_color,
                linewidth=line_width, solid_capstyle="butt", zorder=20,
            ),
            Line2D(
                [line_start, line_start],
                [line_y - cap_half_height, line_y + cap_half_height],
                transform=fig.transFigure, color=line_color,
                linewidth=line_width, solid_capstyle="butt", zorder=20,
            ),
            Line2D(
                [line_end, line_end],
                [line_y - cap_half_height, line_y + cap_half_height],
                transform=fig.transFigure, color=line_color,
                linewidth=line_width, solid_capstyle="butt", zorder=20,
            ),
        ])
        fig.text(
            line_start - 0.006, line_y, format_scale_length(scale_length),
            ha="right", va="center", fontsize=6.8, color=line_color,
        )
        if assembly_label:
            fig.text(
                plot_right + 0.006, line_y, ellipsize(assembly_label, 18),
                ha="left", va="center", fontsize=6.8, color=line_color,
            )

    def draw_coverage_track(
        self, ax, reads: List[AlignedRead], start: int, end: int
    ) -> int:
        """Draw depth with qualifying SNV allele fractions stacked in base colours."""
        span = max(end - start, 1)
        axes_width_pixels = max(int(ax.get_window_extent().width), 1)
        bin_limit = max(
            int(axes_width_pixels * self.styles["coverage_bins_per_pixel"]), 1
        )
        use_binned_coverage = span > bin_limit
        coverage_label = "molecule coverage" if self.molecule_mode else "coverage"
        variant_width = 1.0
        if use_binned_coverage:
            edges, depth = compute_binned_coverage(reads, start, end, bin_limit)
            if depth:
                ax.stairs(
                    depth, edges, baseline=0, fill=True,
                    color=self.visual_colors["coverage"],
                    alpha=self.styles["coverage_alpha"], linewidth=0,
                )
                bin_width = span / len(depth)
                variant_width = max(bin_width * 0.55, 1.0)
                unit = "molecule coverage" if self.molecule_mode else "coverage"
                coverage_label = f"{unit} · {bin_width:.0f} bp/bin (mean)"
        else:
            depth = compute_coverage(reads, start, end)
            positions = []
            for index in range(len(depth)):
                positions.append(start + index + 0.5)
            ax.bar(
                positions, depth, width=1.0, color=self.visual_colors["coverage"],
                alpha=self.styles["coverage_alpha"], linewidth=0,
            )
        max_depth = max(depth) if depth else 0

        if use_binned_coverage:
            evidence_depth, snv_evidence = compute_sparse_snv_evidence(
                reads, start, end,
                min_baseq=self.min_baseq, min_mapq=self.min_variant_mapq,
            )
        else:
            evidence_depth, snv_evidence = compute_snv_evidence(
                reads, start, end,
                min_baseq=self.min_baseq, min_mapq=self.min_variant_mapq,
            )
        labels = []
        for position, base_counts in snv_evidence.items():
            if use_binned_coverage:
                position_depth = evidence_depth.get(position, 0)
            else:
                position_depth = evidence_depth[position - start]
            if position_depth <= 0:
                continue
            max_depth = max(max_depth, position_depth)
            bottom = 0
            for base in "ACGT":
                allele = base_counts.get(base)
                if allele is None or allele.count / position_depth <= self.coverage_vaf_threshold:
                    continue
                ax.bar(
                    position + 0.5, allele.count, width=variant_width, bottom=bottom,
                    color=self.base_colors[base], linewidth=0, zorder=2,
                )
                bottom += allele.count
                labels.append((position, base, allele, position_depth))

        ax.set_ylim(0, max(max_depth, 1) * 1.15)
        ax.set_yticks([0, max(max_depth, 1)])
        ax.tick_params(left=True, labelleft=True, labelsize=6, colors=self.visual_colors["axis"], length=3)
        ax.spines["left"].set_visible(True)
        ax.spines["left"].set_color(self.visual_colors["axis"])
        ax.spines["left"].set_linewidth(0.8)
        ax.text(
            0.005, 0.98, coverage_label, transform=ax.transAxes, fontsize=7,
            color=self.visual_colors["secondary_text"], va="top",
        )

        base_spacing_px = ax.get_window_extent().width / max(end - start, 1)
        if self.show_variant_counts and base_spacing_px >= 5.5:
            for position, base, allele, position_depth in labels:
                mean_baseq = allele.base_quality_sum / allele.count
                mean_mapq = allele.mapq_sum / allele.count
                label = (
                    f"{base} {allele.count}/{position_depth} "
                    f"{allele.count / position_depth:.0%} "
                    f"F{allele.forward}/R{allele.reverse} "
                    f"BQ{mean_baseq:.0f} MQ{mean_mapq:.0f}"
                )
                ax.text(
                    position + 0.5, 0.97, label,
                    transform=ax.get_xaxis_transform(), rotation=90,
                    ha="right", va="center", fontsize=4.8,
                    color=self.visual_colors["primary_text"], clip_on=True, zorder=4,
                    bbox={
                        "facecolor": self.visual_colors["label_background"],
                        "edgecolor": "none", "alpha": 0.68, "pad": 0.15,
                    },
                )
        return max(max_depth, 1)

    def draw_modification_track(
        self, ax, reads: List[AlignedRead], start: int, end: int,
    ) -> None:
        """Draw per-site modified/canonical fractions decoded from MM/ML tags."""
        evidence = compute_modification_evidence(
            reads, start, end,
            minimum_probability=self.min_mod_probability,
            codes=self.modification_codes,
        )
        ax.set_ylim(0, 1.05)
        ax.axhline(0, color=self.visual_colors["axis"], linewidth=0.55, zorder=1)
        for label in sorted(evidence):
            sites = evidence[label]
            positions = sorted(
                position for position, item in sites.items() if item.canonical_depth
            )
            if not positions:
                continue
            fractions = [sites[position].fraction for position in positions]
            color = modification_color(
                label, self.modification_colors, self.chromosome_palette
            )
            self.modification_labels_seen.add(label)
            centers = [position + 0.5 for position in positions]
            ax.vlines(
                centers, 0, fractions, color=color, linewidth=0.7,
                alpha=self.styles["modification_stem_alpha"], zorder=2,
            )
            ax.scatter(
                centers, fractions, s=self.styles["modification_marker_size"],
                marker="o", facecolors=color,
                edgecolors=self.visual_colors["contrast_edge"], linewidths=0.25,
                alpha=self.styles["modification_track_alpha"], zorder=3,
            )
        ax.set_yticks([0, 0.5, 1.0])
        ax.set_yticklabels(["0", ".5", "1"], fontsize=6)
        ax.tick_params(
            left=True, labelleft=True, labelsize=6,
            colors=self.visual_colors["axis"], length=2,
        )
        ax.spines["left"].set_visible(True)
        ax.spines["left"].set_color(self.visual_colors["axis"])
        ax.text(
            0.005, 0.97, "MM/ML modified fraction",
            transform=ax.transAxes, ha="left", va="top",
            fontsize=6.5, color=self.visual_colors["secondary_text"],
            clip_on=True,
        )

    def draw_center_guide(self, ax, start: int, end: int) -> None:
        """Draw the optional IGV-like guide at the exact window midpoint."""
        if not self.show_center_guide:
            return
        ax.axvline(
            start + (end - start) / 2,
            color=self.visual_colors["center_guide"],
            alpha=self.styles["center_guide_alpha"],
            linewidth=self.styles["center_guide_width"],
            linestyle=self.styles["center_guide_line_style"],
            zorder=20,
        )

    def draw_sashimi_track(
        self, ax, reads: List[AlignedRead], start: int, end: int,
        chrom: str = "", reference: Optional[ReferenceWindow] = None,
        genomic_tracks: Optional[List[LoadedAnnotationTrack]] = None,
    ) -> None:
        """Draw classified splice junctions and clustered fusion evidence."""
        junctions = collect_splice_junctions(
            reads, strand_mode=self.sashimi_strand,
            strandness=self.rna_strandness,
            minimum_anchor=self.min_junction_anchor,
        ) if self.show_sashimi else {}
        known = annotated_junctions(genomic_tracks)
        annotation_available = any(
            track.kind in ("bed", "gff", "gff3", "gtf") and track.items
            for track in genomic_tracks or []
        )
        visible = []
        for evidence in junctions.values():
            if (
                evidence.count >= self.min_junction_reads
                and start <= evidence.donor < evidence.acceptor <= end
            ):
                visible.append(evidence)

        fusions = []
        if self.show_fusions:
            fusions = [
                evidence for evidence in collect_fusion_evidence(
                    reads, chrom,
                    breakpoint_tolerance=self.fusion_breakpoint_tolerance,
                    minimum_distance=self.fusion_min_distance,
                    minimum_mapq=self.min_fusion_mapq,
                )
                if evidence.support >= self.min_fusion_reads
                and start <= evidence.local_breakpoint <= end
            ][:6]

        def junction_span_key(item):
            return (item.acceptor - item.donor, item.donor)

        visible.sort(key=junction_span_key, reverse=True)

        if self.sashimi_strand == "split":
            ax.set_ylim(-1.05, 1.28 if self.show_fusions else 1.05)
        else:
            ax.set_ylim(-0.08, 1.28 if self.show_fusions else 1.05)
        ax.axhline(
            0, color=self.visual_colors["axis"], linewidth=0.55, zorder=1
        )
        ax.text(
            -0.012, 0.5,
            "junctions / fusions" if self.show_fusions and self.show_sashimi else
            "fusion evidence" if self.show_fusions else "splice junctions",
            transform=ax.transAxes,
            ha="right", va="center", fontsize=7,
            color=self.visual_colors["sashimi_combined"],
            fontweight="bold", clip_on=False,
        )
        if not visible and not fusions:
            evidence_name = "RNA evidence" if self.show_fusions else "junctions"
            ax.text(
                0.01, 0.60, f"No {evidence_name} above support thresholds",
                transform=ax.transAxes, ha="left", va="center", fontsize=6.5,
                color=self.visual_colors["secondary_text"],
            )
            return

        maximum_span = max(
            (evidence.acceptor - evidence.donor for evidence in visible), default=1
        )
        maximum_count = max((evidence.count for evidence in visible), default=1)
        for evidence in visible:
            donor, acceptor, strand = (
                evidence.donor, evidence.acceptor, evidence.strand
            )
            count = evidence.count
            span = acceptor - donor
            direction = -1 if self.sashimi_strand == "split" and strand == "-" else 1
            height = self.styles["sashimi_arc_height"] * sqrt(span / maximum_span)
            height *= direction
            is_known = is_annotated_junction(chrom, evidence, known)
            status = "known" if is_known else "novel" if annotation_available else "unclassified"
            motif = splice_motif(evidence, reference)
            noncanonical = motif is not None and motif not in CANONICAL_SPLICE_MOTIFS
            if noncanonical:
                color = self.visual_colors["junction_noncanonical"]
                linestyle = ":"
            elif status == "known":
                color = self.visual_colors["junction_annotated"]
                linestyle = "-"
            elif status == "novel":
                color = self.visual_colors["junction_novel"]
                linestyle = "--"
            else:
                color_key = (
                    "sashimi_minus" if direction < 0 else
                    "sashimi_plus" if self.sashimi_strand == "split" else
                    "sashimi_combined"
                )
                color = self.visual_colors[color_key]
                linestyle = "-"
            width_fraction = sqrt(count / maximum_count)
            line_width = self.styles["sashimi_min_line_width"] + width_fraction * (
                self.styles["sashimi_max_line_width"]
                - self.styles["sashimi_min_line_width"]
            )
            vertices = [
                (donor, 0),
                (donor + span * 0.28, height),
                (acceptor - span * 0.28, height),
                (acceptor, 0),
            ]
            path = Path(vertices, [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4])
            ax.add_patch(PathPatch(
                path, facecolor="none", edgecolor=color, linestyle=linestyle,
                linewidth=line_width, alpha=self.styles["sashimi_arc_alpha"],
                capstyle="round", zorder=3,
            ))
            strand_label = strand if self.sashimi_strand == "split" else ""
            status_label = {"known": "K", "novel": "N", "unclassified": "?"}[status]
            label = f"{strand_label}{count}"
            if self.junction_labels in ("status", "full"):
                label += f" {status_label}"
            if self.junction_labels == "full" and motif:
                label += f" {motif}"
            label_y = height + direction * 0.06
            ax.text(
                donor + span / 2, label_y, label,
                ha="center", va="bottom" if direction > 0 else "top",
                fontsize=self.styles["sashimi_label_size"],
                color=color, fontweight="bold", zorder=4,
            )

        fusion_color = self.visual_colors["fusion_split"]
        axis_span = max(end - start, 1)
        for fusion_index, fusion in enumerate(fusions):
            local = fusion.local_breakpoint
            partner_visible = (
                fusion.partner_chrom == chrom
                and start <= fusion.partner_breakpoint <= end
            )
            if partner_visible:
                target = fusion.partner_breakpoint
            else:
                target = end - axis_span * 0.01 if local <= (start + end) / 2 else start + axis_span * 0.01
            direction = 1
            height = max(0.48, 1.08 - fusion_index * 0.11)
            fusion_span = target - local
            vertices = [
                (local, 0),
                (local + fusion_span * 0.28, height),
                (target - fusion_span * 0.18, height),
                (target, height * 0.12),
            ]
            path = Path(
                vertices, [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4]
            )
            ax.add_patch(PathPatch(
                path, facecolor="none", edgecolor=fusion_color,
                linewidth=self.styles["fusion_line_width"] * sqrt(fusion.support),
                alpha=self.styles["fusion_arc_alpha"], capstyle="round", zorder=6,
            ))
            ax.scatter(
                [local], [0], marker="D", s=self.styles["fusion_marker_size"],
                facecolor=fusion_color, edgecolor=self.visual_colors["contrast_edge"],
                linewidth=0.35, zorder=7,
            )
            local_gene = gene_label_at(genomic_tracks, chrom, local)
            partner_label = (
                gene_label_at(genomic_tracks, chrom, fusion.partner_breakpoint)
                if partner_visible else None
            )
            source = local_gene or f"{chrom}:{local + 1:,}"
            partner = partner_label or f"{fusion.partner_chrom}:{fusion.partner_breakpoint + 1:,}"
            support_label = f"S{fusion.split_reads}/P{fusion.spanning_pairs}"
            label_x = local + fusion_span * 0.58
            ax.text(
                label_x, height + direction * 0.045,
                f"{source} → {partner} · {support_label}",
                ha="center", va="bottom" if direction > 0 else "top",
                fontsize=self.styles["fusion_label_size"], color=fusion_color,
                fontweight="bold", clip_on=True, zorder=8,
            )

    def draw_ideogram(
        self,
        ax,
        chrom: str,
        window_start: int,
        window_end: int,
        contig_length: int,
        cytobands: Optional[List[Cytoband]] = None,
    ) -> None:
        """Draw a UCSC-style cytoband ideogram with the current window in red."""
        ax.set_ylim(0, 1)
        ax.set_xlim(window_start, window_end)
        # Use the genomic x-axis itself, not figure-relative coordinates. This
        # locks both chromosome ends to the exact plot boundaries in single,
        # comparison, and mate layouts even when their margins differ.
        chromosome_transform = ax.get_xaxis_transform()
        bar_x = window_start
        bar_y = 0.28
        bar_width = max(window_end - window_start, 1)
        bar_height = 0.44
        chromosome_vertices = None
        neck_position = None
        middle_y = None
        neck_half_height = None
        p_centromeres = []
        q_centromeres = []
        for band in cytobands or []:
            if band.stain != "acen":
                continue
            if band.name.startswith("p"):
                p_centromeres.append(band)
            elif band.name.startswith("q"):
                q_centromeres.append(band)
        if p_centromeres and q_centromeres:
            p_centromere = p_centromeres[0]
            for band in p_centromeres:
                if band.end > p_centromere.end:
                    p_centromere = band
            q_centromere = q_centromeres[0]
            for band in q_centromeres:
                if band.start < q_centromere.start:
                    q_centromere = band
            p_shoulder = bar_x + p_centromere.start / contig_length * bar_width
            q_shoulder = bar_x + q_centromere.end / contig_length * bar_width
            neck_position = (
                p_centromere.end + q_centromere.start
            ) / 2 / contig_length * bar_width + bar_x
            middle_y = bar_y + bar_height / 2
            neck_half_height = bar_height * 0.22
            chromosome_vertices = [
                (bar_x, bar_y),
                (p_shoulder, bar_y),
                (neck_position, middle_y - neck_half_height),
                (q_shoulder, bar_y),
                (bar_x + bar_width, bar_y),
                (bar_x + bar_width, bar_y + bar_height),
                (q_shoulder, bar_y + bar_height),
                (neck_position, middle_y + neck_half_height),
                (p_shoulder, bar_y + bar_height),
                (bar_x, bar_y + bar_height),
            ]
            chromosome_clip = Polygon(
                chromosome_vertices, closed=True,
                transform=chromosome_transform, facecolor=self.visual_colors["ideogram"],
                edgecolor=self.visual_colors["axis"], linewidth=0.6, zorder=2,
            )
        else:
            chromosome_clip = Rectangle(
                (bar_x, bar_y), bar_width, bar_height,
                transform=chromosome_transform, facecolor=self.visual_colors["ideogram"],
                edgecolor=self.visual_colors["axis"], linewidth=0.6, zorder=2,
            )
        ax.add_patch(chromosome_clip)

        for band in cytobands or []:
            band_start = min(max(band.start, 0), contig_length)
            band_end = min(max(band.end, band_start), contig_length)
            if band_end <= band_start:
                continue
            x0 = bar_x + band_start / contig_length * bar_width
            x1 = bar_x + band_end / contig_length * bar_width
            color = self.cytoband_colors.get(band.stain, self.visual_colors["ideogram"])
            patch = Rectangle(
                (x0, bar_y), x1 - x0, bar_height,
                transform=chromosome_transform, facecolor=color,
                edgecolor=self.visual_colors["cytoband_edge"], linewidth=0.15, zorder=2.1,
            )
            patch.set_clip_path(chromosome_clip)
            ax.add_patch(patch)

        if chromosome_vertices:
            bridge_width = max(bar_width * 0.004, 0.4)
            centromere_bridge = Rectangle(
                (neck_position - bridge_width / 2, middle_y - neck_half_height),
                bridge_width, neck_half_height * 2,
                transform=chromosome_transform,
                facecolor=self.visual_colors["centromere"], edgecolor="none", zorder=2.15,
            )
            centromere_bridge.set_clip_path(chromosome_clip)
            ax.add_patch(centromere_bridge)

        if cytobands:
            if chromosome_vertices:
                outline = Polygon(
                    chromosome_vertices, closed=True,
                    transform=chromosome_transform, facecolor="none",
                    edgecolor=self.visual_colors["axis"], linewidth=0.6, zorder=2.2,
                )
            else:
                outline = Rectangle(
                    (bar_x, bar_y), bar_width, bar_height,
                    transform=chromosome_transform, facecolor="none",
                    edgecolor=self.visual_colors["axis"], linewidth=0.6, zorder=2.2,
                )
            ax.add_patch(outline)

        clamped_start = min(max(window_start, 0), contig_length)
        clamped_end = min(max(window_end, clamped_start), contig_length)
        relative_start = clamped_start / contig_length
        relative_width = max((clamped_end - clamped_start) / contig_length, 0.0)
        # Base-pair windows on chromosome-scale bars would otherwise disappear.
        marker_width = max(relative_width * bar_width, bar_width * 0.004)
        marker_center = bar_x + (relative_start + relative_width / 2) * bar_width
        marker_x = min(max(marker_center - marker_width / 2, bar_x), bar_x + bar_width - marker_width)
        ax.add_patch(Rectangle(
            (marker_x, bar_y - 0.08), marker_width, bar_height + 0.16,
            transform=chromosome_transform, facecolor=self.visual_colors["ideogram_window"],
            edgecolor=self.visual_colors["contrast_edge"], linewidth=0.35, zorder=3,
        ))
        ax.text(
            -0.012, 0.5, chrom, transform=ax.transAxes, ha="right", va="center",
            fontsize=7, color=self.visual_colors["primary_text"], fontweight="bold",
            clip_on=False,
        )
        ax.text(
            1.012, 0.5, f"{contig_length / 1_000_000:.1f} Mb",
            transform=ax.transAxes, ha="left", va="center", fontsize=6.5,
            color=self.visual_colors["secondary_text"], clip_on=False,
        )

    def draw_annotation_track(
        self,
        ax,
        track: Any,
        window_start: int,
        window_end: int,
        shared_row_count: Optional[int] = None,
    ) -> None:
        """Draw a UCSC-like BED or transcript annotation track."""
        if isinstance(track, LoadedPluginTrack):
            self.draw_plugin_track(ax, track, window_start, window_end)
            return
        if track.kind in HIC_TRACK_FORMATS:
            self.draw_hic_track(ax, track, window_start, window_end)
            return
        if track.display_mode == "density":
            self.draw_density_track(ax, track, window_start, window_end)
            return
        if track.kind in PEAK_TRACK_FORMATS or track.kind in SIGNAL_TRACK_FORMATS:
            self.draw_peak_track(ax, track, window_start, window_end)
            return
        if track.kind in CNV_TRACK_FORMATS:
            self.draw_cnv_track(ax, track, window_start, window_end)
            return
        if track.kind in BAF_TRACK_FORMATS:
            self.draw_baf_track(ax, track, window_start, window_end)
            return
        row_count = max(shared_row_count or len(track.rows), 1)
        ax.set_ylim(row_count, 0)
        margin_in = min(
            max(0.70, 0.25 + len(track.label) * 0.065), self.fig_width * 0.22
        )
        label_capacity = max(5, int((margin_in - 0.20) / 0.065))
        ax.text(
            -0.012, 0.5, ellipsize(track.label, label_capacity),
            transform=ax.transAxes, ha="right", va="center",
            fontsize=7, color=track.color, fontweight="bold", clip_on=False,
        )
        if not track.rows:
            ax.text(
                0.01, 0.5, "No features", transform=ax.transAxes, ha="left", va="center",
                fontsize=6.5, color=self.visual_colors["axis"],
            )
            return

        for row_index, row in enumerate(track.rows):
            # Reserve the upper quarter of each row for the feature name.
            center = row_index + 0.64
            for item in row:
                line_start = max(item.start, window_start)
                line_end = min(item.end, window_end)
                if line_start >= line_end:
                    continue
                ax.plot(
                    [line_start, line_end], [center, center], color=track.color,
                    linewidth=(
                        self.styles["primary_gene_line_width"] if item.primary_rank is not None
                        else self.styles["gene_line_width"]
                    ),
                    zorder=2, solid_capstyle="butt",
                )
                for block_start, block_end in item.blocks:
                    lo, hi = max(block_start, window_start), min(block_end, window_end)
                    if lo < hi:
                        ax.add_patch(Rectangle(
                            (lo, center - 0.23), hi - lo, 0.46,
                            facecolor=track.color, edgecolor=track.color,
                            linewidth=self.styles["annotation_edge_width"], zorder=3,
                        ))
                for utr_start, utr_end in item.utrs:
                    lo, hi = max(utr_start, window_start), min(utr_end, window_end)
                    if lo < hi:
                        ax.add_patch(Rectangle(
                            (lo, center - 0.12), hi - lo, 0.24,
                            facecolor=track.color, edgecolor=track.color,
                            linewidth=self.styles["annotation_edge_width"], zorder=3,
                        ))

                # Repeated small arrows on introns make transcript direction
                # readable without competing with exon blocks.
                merged_intervals: List[List[int]] = []
                for feature_start, feature_end in sorted(item.blocks + item.utrs):
                    if feature_end <= feature_start:
                        continue
                    if not merged_intervals or feature_start > merged_intervals[-1][1]:
                        merged_intervals.append([feature_start, feature_end])
                    else:
                        merged_intervals[-1][1] = max(
                            merged_intervals[-1][1], feature_end
                        )
                introns = []
                intron_cursor = item.start
                for feature_start, feature_end in merged_intervals:
                    if intron_cursor < feature_start:
                        introns.append((intron_cursor, feature_start))
                    intron_cursor = max(intron_cursor, feature_end)
                if intron_cursor < item.end:
                    introns.append((intron_cursor, item.end))
                if item.strand in ("+", "-") and introns:
                    axes_width_px = max(ax.get_window_extent().width, 1)
                    arrow_spacing = max(
                        (window_end - window_start) / max(
                            axes_width_px / self.styles["gene_arrow_spacing_px"], 1
                        ),
                        1.0,
                    )
                    marker = ">" if item.strand == "+" else "<"
                    for intron_start, intron_end in introns:
                        lo = max(intron_start, window_start)
                        hi = min(intron_end, window_end)
                        if lo >= hi:
                            continue
                        visible_intron_length = hi - lo
                        arrow_count = max(
                            1, int(visible_intron_length / arrow_spacing)
                        )
                        position_step = visible_intron_length / (arrow_count + 1)
                        for arrow_index in range(arrow_count):
                            position = lo + position_step * (arrow_index + 1)
                            ax.plot(
                                position, center, marker=marker,
                                markersize=self.styles["gene_arrow_size"],
                                color=track.color, markeredgewidth=0, zorder=4,
                            )
                elif item.strand == "+" and item.end <= window_end:
                    ax.plot(item.end, center, marker=">",
                            markersize=self.styles["gene_arrow_size"],
                            color=track.color, markeredgewidth=0, zorder=4)
                elif item.strand == "-" and item.start >= window_start:
                    ax.plot(item.start, center, marker="<",
                            markersize=self.styles["gene_arrow_size"],
                            color=track.color, markeredgewidth=0, zorder=4)

                visible_fraction = (line_end - line_start) / max(window_end - window_start, 1)
                name_capacity = int(visible_fraction * 105)
                if track.display_mode == "collapse":
                    item_label = item.group_label or item.name
                else:
                    gene_label = item.group_label or item.group
                    transcript_label = item.transcript_label or item.name
                    if (
                        gene_label and transcript_label
                        and gene_label != transcript_label
                    ):
                        item_label = f"{gene_label} · {transcript_label}"
                    else:
                        item_label = gene_label or transcript_label
                if item.primary_rank is not None:
                    item_label += " ★"
                display_name = ellipsize(item_label, name_capacity) if name_capacity >= 4 else ""
                if display_name:
                    ax.text(
                        line_start, row_index + 0.05, display_name, ha="left", va="top",
                        fontsize=5.5, color=track.color, clip_on=True, zorder=5,
                    )

    def draw_hic_track(
        self, ax, track: LoadedAnnotationTrack, window_start: int, window_end: int
    ) -> None:
        """Draw called TAD domains or BEDPE Hi-C contacts in genomic space."""
        ax.set_ylim(0, 1.05)
        ax.set_yticks([])
        ax.text(
            -0.012, 0.5, ellipsize(track.label, 28),
            transform=ax.transAxes, ha="right", va="center",
            fontsize=7, color=track.color, fontweight="bold", clip_on=False,
        )
        if not track.items:
            ax.text(
                0.01, 0.5, "No domains" if track.kind in TAD_TRACK_FORMATS else "No contacts",
                transform=ax.transAxes, ha="left", va="center",
                fontsize=6.5, color=self.visual_colors["axis"],
            )
            return

        scored = [item.value for item in track.items if item.value is not None]
        score_min = min(scored) if scored else 0.0
        score_max = max(scored) if scored else 1.0

        def score_fraction(item: AnnotationItem) -> float:
            if item.value is None or score_max <= score_min:
                return 1.0
            return max(0.0, min(1.0, (item.value - score_min) / (score_max - score_min)))

        if (
            track.kind in HIC_LOOP_TRACK_FORMATS
            and track.display_mode == "triangle"
        ):
            self.draw_hic_contact_map(
                ax, track, window_start, window_end, score_fraction
            )
            return

        if track.kind in TAD_TRACK_FORMATS:
            for item in sorted(track.items, key=lambda value: (value.start, value.end)):
                left = max(item.start, window_start)
                right = min(item.end, window_end)
                if right <= left:
                    continue
                midpoint = (item.start + item.end) / 2
                visible_midpoint = max(left, min(right, midpoint))
                height = 0.68 + 0.24 * score_fraction(item)
                triangle = Polygon(
                    [(left, 0.08), (visible_midpoint, height), (right, 0.08)],
                    closed=True, facecolor=track.color,
                    edgecolor=track.color,
                    linewidth=self.styles["tad_line_width"],
                    alpha=self.styles["tad_fill_alpha"], zorder=2,
                )
                ax.add_patch(triangle)
                for boundary in (item.start, item.end):
                    if window_start <= boundary <= window_end:
                        ax.plot(
                            [boundary, boundary], [0.04, height], color=track.color,
                            linewidth=self.styles["tad_line_width"],
                            alpha=self.styles["tad_boundary_alpha"], zorder=3,
                        )
                if item.name and right - left >= (window_end - window_start) * 0.04:
                    ax.text(
                        (left + right) / 2, min(height + 0.025, 1.0),
                        ellipsize(item.name, 28),
                        ha="center", va="bottom", fontsize=5.5,
                        color=track.color, clip_on=True, zorder=4,
                    )
            ax.text(
                0.005, 0.98, "TAD domains · boundary guides",
                transform=ax.transAxes, ha="left", va="top", fontsize=5.5,
                color=self.visual_colors["secondary_text"],
            )
            return

        window_span = max(window_end - window_start, 1)
        for item in track.items:
            if item.start2 is None or item.end2 is None:
                continue
            center1 = (item.start + item.end) / 2
            center2 = (item.start2 + item.end2) / 2
            line_width = self.styles["hic_loop_min_width"] + score_fraction(item) * (
                self.styles["hic_loop_max_width"] - self.styles["hic_loop_min_width"]
            )
            if item.group == item.chrom2:
                left_center, right_center = sorted((center1, center2))
                visible_left = max(window_start, left_center)
                visible_right = min(window_end, right_center)
                if visible_right <= visible_left:
                    continue
                distance_fraction = min(1.0, abs(center2 - center1) / window_span)
                arc_height = 0.16 + 0.76 * sqrt(distance_fraction)
                path = Path(
                    [(visible_left, 0.10),
                     ((visible_left + visible_right) / 2, arc_height),
                     (visible_right, 0.10)],
                    [Path.MOVETO, Path.CURVE3, Path.CURVE3],
                )
                ax.add_patch(PathPatch(
                    path, fill=False, edgecolor=track.color,
                    linewidth=line_width, alpha=self.styles["hic_loop_alpha"],
                    capstyle="round", zorder=2,
                ))
                for anchor_start, anchor_end in (
                    (item.start, item.end), (item.start2, item.end2)
                ):
                    left = max(anchor_start, window_start)
                    right = min(anchor_end, window_end)
                    if right > left:
                        ax.add_patch(Rectangle(
                            (left, 0.055), right - left, 0.09,
                            facecolor=track.color, edgecolor="none",
                            alpha=self.styles["hic_anchor_alpha"], zorder=3,
                        ))
                if item.name and visible_right - visible_left >= window_span * 0.04:
                    ax.text(
                        (visible_left + visible_right) / 2, min(arc_height + 0.035, 1.0),
                        ellipsize(item.name, 24), ha="center", va="bottom",
                        fontsize=5.2, color=track.color, clip_on=True,
                    )
                continue

            visible_anchors = []
            if (
                item.group == track.chrom
                and item.end > window_start and item.start < window_end
            ):
                visible_anchors.append((item.start, item.end, item.chrom2))
            if (
                item.chrom2 == track.chrom
                and item.end2 > window_start and item.start2 < window_end
            ):
                visible_anchors.append((item.start2, item.end2, item.group))
            for anchor_start, anchor_end, partner_chrom in visible_anchors:
                left = max(anchor_start, window_start)
                right = min(anchor_end, window_end)
                center = (left + right) / 2
                ax.add_patch(Rectangle(
                    (left, 0.055), right - left, 0.09,
                    facecolor=track.color, edgecolor="none",
                    alpha=self.styles["hic_anchor_alpha"], zorder=3,
                ))
                ax.plot(
                    center, 0.22, marker="^", markersize=4.2,
                    color=track.color, markeredgewidth=0, zorder=3,
                )
                ax.text(
                    center, 0.28, f"to {partner_chrom}", ha="center", va="bottom",
                    fontsize=5.2, color=track.color, clip_on=True,
                )
        ax.text(
            0.005, 0.98, "Hi-C contacts · BEDPE anchors",
            transform=ax.transAxes, ha="left", va="top", fontsize=5.5,
            color=self.visual_colors["secondary_text"],
        )

    def draw_hic_contact_map(
        self, ax, track: LoadedAnnotationTrack,
        window_start: int, window_end: int, score_fraction,
    ) -> None:
        """Render scored cis BEDPE bins as a rotated triangular contact map."""
        cis_items = [
            item for item in track.items
            if item.group == item.chrom2 == track.chrom
            and item.start2 is not None and item.end2 is not None
        ]
        if not cis_items:
            ax.text(
                0.01, 0.5, "No cis contacts", transform=ax.transAxes,
                ha="left", va="center", fontsize=6.5,
                color=self.visual_colors["axis"],
            )
            return

        low_rgb = to_rgb(self.track_colors["hic_contact_low"])
        high_rgb = to_rgb(track.color)
        gamma = self.styles["hic_contact_gamma"]
        window_span = max(window_end - window_start, 1)

        def contact_color(item: AnnotationItem):
            intensity = score_fraction(item) ** gamma
            return tuple(
                low + (high - low) * intensity
                for low, high in zip(low_rgb, high_rgb)
            )

        for item in sorted(cis_items, key=score_fraction):
            anchor1 = (item.start, item.end)
            anchor2 = (item.start2, item.end2)
            if sum(anchor1) > sum(anchor2):
                anchor1, anchor2 = anchor2, anchor1
            start1, end1 = anchor1
            start2, end2 = anchor2
            matrix_corners = (
                (start1, start2), (end1, start2),
                (end1, end2), (start1, end2),
            )
            polygon_points = [
                (
                    (position1 + position2) / 2,
                    0.05 + 0.90 * abs(position2 - position1) / window_span,
                )
                for position1, position2 in matrix_corners
            ]
            projected_x = [point[0] for point in polygon_points]
            if max(projected_x) < window_start or min(projected_x) > window_end:
                continue
            color = contact_color(item)
            ax.add_patch(Polygon(
                polygon_points, closed=True, facecolor=color, edgecolor=color,
                linewidth=self.styles["hic_contact_cell_edge_width"],
                alpha=self.styles["hic_contact_map_alpha"], zorder=2,
            ))

        ax.plot(
            [window_start, window_end], [0.05, 0.05],
            color=track.color, linewidth=0.45, alpha=0.45, zorder=3,
        )
        ax.text(
            0.005, 0.98, "Triangular Hi-C contact map · score intensity",
            transform=ax.transAxes, ha="left", va="top", fontsize=5.5,
            color=self.visual_colors["secondary_text"],
        )

    def draw_plugin_track(
        self, ax, track: LoadedPluginTrack,
        window_start: int, window_end: int,
    ) -> None:
        """Render an external track through the public versioned canvas API."""
        ax.set_ylim(0, 1)
        ax.text(
            0.005, 0.97, ellipsize(track.label, max(18, int(self.fig_width * 12))),
            transform=ax.transAxes, ha="left", va="top",
            fontsize=7, color=track.color, fontweight="bold", clip_on=True,
        )
        canvas = TrackCanvas(
            ax, track.region, track.color,
            colors=self.visual_colors,
        )
        try:
            track.plugin.render(
                canvas, track.payload, track.region, track.options
            )
        except Exception as exc:
            raise TrackPluginError(
                f"Track plugin '{track.plugin.name}' failed while rendering "
                f"{track.region.chrom}:{window_start + 1}-{window_end}: {exc}"
            ) from exc

    def draw_peak_track(
        self, ax, track: LoadedAnnotationTrack, window_start: int, window_end: int
    ) -> None:
        """Draw called peaks as blocks or quantitative signal as one profile."""
        continuous_signal = track.kind in SIGNAL_TRACK_FORMATS
        largest_value = 1.0
        for item in track.items:
            if item.value is not None and item.value > largest_value:
                largest_value = item.value
        configured_maximum = self.styles["signal_y_max"] if continuous_signal else 0
        display_maximum = (
            configured_maximum if configured_maximum > 0
            else max(largest_value, 1.0)
        )
        value_limit = display_maximum * 1.12
        ax.set_ylim(0, value_limit)
        ax.axhline(0, color=self.visual_colors["axis"], linewidth=0.55, zorder=1)

        if continuous_signal:
            positions = [window_start]
            values = [0.0]
            cursor = window_start
            ordered_items = sorted(
                track.items, key=lambda item: (item.start, item.end)
            )
            for item in ordered_items:
                lo = max(item.start, window_start, cursor)
                hi = min(item.end, window_end)
                if lo >= hi:
                    continue
                if lo > cursor:
                    positions.extend([cursor, lo])
                    values.extend([0.0, 0.0])
                value = item.value if item.value is not None else 0.0
                positions.extend([lo, hi])
                values.extend([value, value])
                cursor = hi
            positions.extend([cursor, window_end])
            values.extend([0.0, 0.0])
            if len(positions) > 2:
                ax.fill_between(
                    positions, 0, values, facecolor=track.color,
                    alpha=self.styles["signal_fill_alpha"], linewidth=0,
                    zorder=2,
                )
                ax.plot(
                    positions, values, color=track.color,
                    linewidth=self.styles["signal_line_width"], zorder=3,
                    solid_capstyle="butt", solid_joinstyle="miter",
                )
        else:
            for item in track.items:
                lo, hi = max(item.start, window_start), min(item.end, window_end)
                if lo >= hi:
                    continue
                value = item.value if item.value is not None else 1.0
                ax.add_patch(Rectangle(
                    (lo, 0), hi - lo, value,
                    facecolor=track.color, edgecolor=track.color,
                    linewidth=0.35, alpha=self.styles["peak_fill_alpha"], zorder=2,
                ))
                if item.summit is not None and window_start <= item.summit < window_end:
                    ax.plot(
                        [item.summit, item.summit], [0, value], color=track.color,
                        linewidth=self.styles["peak_summit_width"], zorder=3,
                    )
                    ax.plot(
                        item.summit, value, marker="v", markersize=3.0,
                        color=track.color, markeredgewidth=0, zorder=4,
                    )

        margin_in = min(
            max(0.70, 0.25 + len(track.label) * 0.065), self.fig_width * 0.22
        )
        label_capacity = max(5, int((margin_in - 0.20) / 0.065))
        ax.text(
            -0.012, 0.5, ellipsize(track.label, label_capacity),
            transform=ax.transAxes, ha="right", va="center", fontsize=7,
            color=track.color, fontweight="bold", clip_on=False,
        )
        if not track.items:
            ax.text(
                0.01, 0.5, "No peaks / signal", transform=ax.transAxes,
                ha="left", va="center", fontsize=6.5,
                color=self.visual_colors["axis"],
            )
        ax.set_yticks([0, display_maximum])
        ax.tick_params(
            right=True, labelright=True, left=False, labelleft=False,
            labelsize=6, colors=self.visual_colors["axis"], length=2,
        )
        ax.spines["right"].set_visible(True)
        ax.spines["right"].set_color(self.visual_colors["axis"])
        ax.spines["right"].set_linewidth(0.6)
        ax.yaxis.set_label_position("right")
        ax.set_ylabel(
            "normalized signal" if continuous_signal else "signalValue",
            rotation=90, labelpad=18, fontsize=5.5,
            color=self.visual_colors["secondary_text"], va="center",
        )

    def draw_density_track(
        self, ax, track: LoadedAnnotationTrack, window_start: int, window_end: int
    ) -> None:
        """Draw compact binned interval-overlap density for any annotation source."""
        axes_width = max(int(ax.get_window_extent().width), 20)
        bin_count = min(max(axes_width // 3, 20), 400)
        positions, densities = compute_feature_density(
            track.items, window_start, window_end, bin_count
        )
        largest_density = max(densities, default=0)
        display_maximum = max(largest_density, 1)
        value_limit = display_maximum * 1.12
        ax.set_ylim(0, value_limit)
        if positions:
            ax.fill_between(
                positions, densities, 0, step="mid", color=track.color,
                alpha=self.styles["density_fill_alpha"], linewidth=0, zorder=2,
            )
            ax.plot(
                positions, densities, drawstyle="steps-mid", color=track.color,
                linewidth=0.75, zorder=3,
            )
        ax.axhline(0, color=self.visual_colors["axis"], linewidth=0.55, zorder=1)

        margin_in = min(
            max(0.70, 0.25 + len(track.label) * 0.065), self.fig_width * 0.22
        )
        label_capacity = max(5, int((margin_in - 0.20) / 0.065))
        ax.text(
            -0.012, 0.5, ellipsize(track.label, label_capacity),
            transform=ax.transAxes, ha="right", va="center", fontsize=7,
            color=track.color, fontweight="bold", clip_on=False,
        )
        ax.text(
            0.005, 0.95, "density", transform=ax.transAxes,
            ha="left", va="top", fontsize=5.5,
            color=self.visual_colors["secondary_text"], clip_on=True,
        )
        if not track.items:
            ax.text(
                0.01, 0.5, "No features", transform=ax.transAxes,
                ha="left", va="center", fontsize=6.5,
                color=self.visual_colors["axis"],
            )
        ax.set_yticks([0, display_maximum])
        ax.tick_params(
            right=True, labelright=True, left=False, labelleft=False,
            labelsize=6, colors=self.visual_colors["axis"], length=2,
        )
        ax.spines["right"].set_visible(True)
        ax.spines["right"].set_color(self.visual_colors["axis"])
        ax.spines["right"].set_linewidth(0.6)
        ax.yaxis.set_label_position("right")
        ax.set_ylabel(
            "features/bin", rotation=90, labelpad=18, fontsize=5.5,
            color=self.visual_colors["secondary_text"], va="center",
        )

    def draw_cnv_track(
        self, ax, track: LoadedAnnotationTrack, window_start: int, window_end: int
    ) -> None:
        """Draw segmented or binned log2 copy-number values around a zero baseline."""
        largest_value = 0.0
        for item in track.items:
            if item.value is not None and abs(item.value) > largest_value:
                largest_value = abs(item.value)
        value_limit = max(0.5, ceil(largest_value * 1.15 * 2) / 2)
        ax.set_ylim(-value_limit, value_limit)
        ax.axhline(0, color=self.visual_colors["axis"], linewidth=0.65, zorder=1)

        use_sign_colors = track.color_by_sign
        for item in track.items:
            if item.value is None:
                continue
            lo, hi = max(item.start, window_start), min(item.end, window_end)
            if lo >= hi:
                continue
            if use_sign_colors:
                color = (
                    self.visual_colors["cnv_gain"] if item.value > 0
                    else self.visual_colors["cnv_loss"] if item.value < 0
                    else track.color
                )
            else:
                color = track.color
            ax.add_patch(Rectangle(
                (lo, min(0, item.value)), hi - lo, abs(item.value),
                facecolor=color, edgecolor="none",
                alpha=self.styles["cnv_fill_alpha"], zorder=2,
            ))
            ax.plot(
                [lo, hi], [item.value, item.value], color=color,
                linewidth=1.25, solid_capstyle="butt", zorder=3,
            )

        margin_in = min(
            max(0.70, 0.25 + len(track.label) * 0.065), self.fig_width * 0.22
        )
        label_capacity = max(5, int((margin_in - 0.20) / 0.065))
        ax.text(
            -0.012, 0.5, ellipsize(track.label, label_capacity),
            transform=ax.transAxes, ha="right", va="center", fontsize=7,
            color=track.color, fontweight="bold", clip_on=False,
        )
        sample_names = []
        seen_samples = set()
        for item in track.items:
            if item.sample and item.sample not in seen_samples:
                seen_samples.add(item.sample)
                sample_names.append(item.sample)
        if sample_names:
            sample_label = ", ".join(sample_names[:3])
            if len(sample_names) > 3:
                sample_label += f" +{len(sample_names) - 3}"
            ax.text(
                0.005, 0.97, ellipsize(sample_label, max(12, int(self.fig_width * 8))),
                transform=ax.transAxes, ha="left", va="top", fontsize=6,
                color=self.visual_colors["secondary_text"], clip_on=True,
            )
        elif not track.items:
            ax.text(
                0.01, 0.5, "No CNV data", transform=ax.transAxes,
                ha="left", va="center", fontsize=6.5, color=self.visual_colors["axis"],
            )
        ax.set_yticks([-value_limit, 0, value_limit])
        ax.set_yticklabels([
            f"{-value_limit:g}", "0", f"{value_limit:g}",
        ])
        ax.tick_params(
            right=True, labelright=True, left=False, labelleft=False,
            labelsize=6, colors=self.visual_colors["axis"], length=2,
        )
        ax.spines["right"].set_visible(True)
        ax.spines["right"].set_color(self.visual_colors["axis"])
        ax.spines["right"].set_linewidth(0.6)
        ax.yaxis.set_label_position("right")
        ax.set_ylabel(
            "log2", rotation=90, labelpad=18, fontsize=5.5,
            color=self.visual_colors["secondary_text"], va="center",
        )

    def draw_baf_track(
        self, ax, track: LoadedAnnotationTrack, window_start: int, window_end: int
    ) -> None:
        """Draw heterozygous-SNV B-allele fractions on a zero-to-one scale."""
        ax.set_ylim(0, 1)
        ax.axhline(0.5, color=self.visual_colors["axis"], linewidth=0.65, linestyle="--", zorder=1)
        positions = []
        values = []
        for item in track.items:
            if item.value is not None:
                positions.append(item.start + 0.5)
                values.append(item.value)
        if positions:
            ax.scatter(
                positions, values, s=11, color=track.color,
                edgecolors=self.visual_colors["contrast_edge"], linewidths=0.25,
                alpha=self.styles["baf_alpha"], zorder=3,
            )

        margin_in = min(
            max(0.70, 0.25 + len(track.label) * 0.065), self.fig_width * 0.22
        )
        label_capacity = max(5, int((margin_in - 0.20) / 0.065))
        ax.text(
            -0.012, 0.5, ellipsize(track.label, label_capacity),
            transform=ax.transAxes, ha="right", va="center", fontsize=7,
            color=track.color, fontweight="bold", clip_on=False,
        )
        sample_names = []
        seen_samples = set()
        for item in track.items:
            if item.sample and item.sample not in seen_samples:
                seen_samples.add(item.sample)
                sample_names.append(item.sample)
        if sample_names:
            ax.text(
                0.005, 0.97, ellipsize(", ".join(sample_names), max(12, int(self.fig_width * 8))),
                transform=ax.transAxes, ha="left", va="top", fontsize=6,
                color=self.visual_colors["secondary_text"], clip_on=True,
            )
        elif not track.items:
            ax.text(
                0.01, 0.5, "No heterozygous SNPs with AD/AF",
                transform=ax.transAxes, ha="left", va="center",
                fontsize=6.5, color=self.visual_colors["axis"],
            )
        ax.set_yticks([0, 0.5, 1])
        ax.tick_params(
            right=True, labelright=True, left=False, labelleft=False,
            labelsize=6, colors=self.visual_colors["axis"], length=2,
        )
        ax.spines["right"].set_visible(True)
        ax.spines["right"].set_color(self.visual_colors["axis"])
        ax.spines["right"].set_linewidth(0.6)
        ax.yaxis.set_label_position("right")
        ax.set_ylabel(
            "BAF", rotation=90, labelpad=18, fontsize=5.5,
            color=self.visual_colors["secondary_text"], va="center",
        )

    def draw_reference_track(
        self,
        ax,
        reference: ReferenceWindow,
        window_start: int,
        window_end: int,
        available_width_in: float,
    ) -> None:
        """Draw one lightly coloured cell per FASTA base, with letters when legible."""
        ax.set_ylim(0, 1)
        span = max(window_end - window_start, 1)
        show_letters = available_width_in * 72 / span >= 6.2
        for pos in range(window_start, window_end):
            base = reference.base_at(pos) or "N"
            color = self.base_colors.get(base, self.base_colors["N"])
            ax.add_patch(Rectangle(
                (pos, 0.08), 1, 0.84, facecolor=color,
                alpha=self.styles["reference_base_alpha"],
                edgecolor=self.visual_colors["contrast_edge"], linewidth=0.25, zorder=2,
            ))
            if show_letters:
                ax.text(
                    pos + 0.5, 0.5, base, ha="center", va="center",
                    fontsize=7, color=color, fontweight="bold", zorder=3,
                    clip_on=True,
                )

    def draw_legends(
        self, fig, fig_height: float, plot_left: float = 0.05,
        plot_right: float = 0.95,
    ) -> list:
        """Draw responsive topic cards below the genomic plot."""
        if self.molecule_mode:
            alignment_handles = [
                Patch(facecolor=self.molecule_colors["singleton"], edgecolor="none", label="Singleton"),
                Patch(facecolor=self.molecule_colors["consensus"], edgecolor="none", label="Consensus family"),
                Patch(facecolor=self.molecule_colors["duplex"], edgecolor="none", label="Duplex family"),
            ]
        elif self.long_read_coloring:
            alignment_handles = [
                Patch(facecolor=self.long_read_colors["forward"], edgecolor="none", label="Forward"),
                Patch(facecolor=self.long_read_colors["reverse"], edgecolor="none", label="Reverse"),
                Patch(facecolor=self.long_read_colors["supplementary"], edgecolor="none", label="Supplementary"),
            ]
        else:
            alignment_handles = [
                Patch(facecolor=self.alignment_colors["normal"], edgecolor="none", label="Normal / concordant"),
            ]
        if self.view_as_pairs:
            alignment_handles.append(
                Line2D([0], [0], color=self.visual_colors["secondary_text"], lw=1.0, label="Mate link")
            )
        event_handles = [
            Patch(facecolor=self.visual_colors["insertion"], edgecolor="none", label="Insertion"),
            Line2D([0], [0], color=self.visual_colors["deletion"], lw=1.5, label="Deletion"),
        ]
        insert_size_handles = []
        pair_geometry_handles = []
        if (
            self.pair_colors
            and self.haplotype_view == "none"
            and self.tag_view == "none"
        ):
            insert_size_handles = [
                Patch(facecolor=self.alignment_colors["large_insert"], edgecolor="none", label=PAIR_CATEGORY_LABELS["large_insert"]),
                Patch(facecolor=self.alignment_colors["small_insert"], edgecolor="none", label=PAIR_CATEGORY_LABELS["small_insert"]),
            ]
            pair_geometry_handles = [
                Patch(facecolor=self.alignment_colors["ff"], edgecolor="none", label=PAIR_CATEGORY_LABELS["ff"]),
                Patch(facecolor=self.alignment_colors["rr"], edgecolor="none", label=PAIR_CATEGORY_LABELS["rr"]),
                Patch(facecolor=self.alignment_colors["everted"], edgecolor="none", label=PAIR_CATEGORY_LABELS["everted"]),
            ]
            if self.alignment_colors["interchrom"]:
                pair_geometry_handles.append(Patch(
                    facecolor=self.alignment_colors["interchrom"], edgecolor="none",
                    label=PAIR_CATEGORY_LABELS["interchrom"],
                ))
            elif self.interchrom_mate_colors:
                mate_chromosomes = sorted(self.interchrom_mate_colors.items())
                if len(mate_chromosomes) <= MAX_EXPLICIT_MATE_CHROMOSOMES:
                    for mate_chrom, mate_color in mate_chromosomes:
                        pair_geometry_handles.append(Patch(
                            facecolor=mate_color, edgecolor="none",
                            label=f"Mate {mate_chrom}",
                        ))
                else:
                    # Preserve mate-chromosome colouring in the reads, but keep
                    # rearrangement-heavy legends to one bounded summary item.
                    pair_geometry_handles.append(Patch(
                        facecolor=mate_chromosomes[0][1], edgecolor="none",
                        label=(
                            f"Inter-chromosomal ({len(mate_chromosomes)} chromosomes)"
                        ),
                    ))
            else:
                pair_geometry_handles.append(Patch(
                    facecolor=chrom_color(
                        "chr1", self.chromosome_palette,
                        colors=self.chromosome_colors,
                    ),
                    edgecolor="none", label="Inter-chromosomal (mate colour)",
                ))
            alignment_handles.extend(pair_geometry_handles)
        haplotype_handles = []
        if self.haplotype_view in ("color", "split"):
            haplotype_handles = [
                Patch(facecolor=haplotype_color("1", self.haplotype_colors, self.chromosome_palette), edgecolor="none", label="HP 1"),
                Patch(facecolor=haplotype_color("2", self.haplotype_colors, self.chromosome_palette), edgecolor="none", label="HP 2"),
                Patch(facecolor=haplotype_color("3", self.haplotype_colors, self.chromosome_palette), edgecolor="none", label="Other HP"),
                Patch(facecolor=haplotype_color(None, self.haplotype_colors, self.chromosome_palette), edgecolor="none", label="Untagged"),
            ]
        tag_handles = []
        if self.tag_view in ("color", "split") and self.tag_value_colors:
            labels = sorted(
                self.tag_value_colors,
                key=lambda label: (
                    label == "untagged",
                    0 if label.isdigit() else 1,
                    int(label) if label.isdigit() else label,
                ),
            )
            omitted_count = 0
            if len(labels) > MAX_EXPLICIT_TAG_VALUES:
                labels = sorted(
                    labels,
                    key=lambda label: (-self.tag_value_counts.get(label, 0), label),
                )[:MAX_EXPLICIT_TAG_VALUES - 1]
                omitted_count = len(self.tag_value_colors) - len(labels)
            for label in labels:
                display_label = "Untagged" if label == "untagged" else label
                tag_handles.append(Patch(
                    facecolor=self.tag_value_colors[label], edgecolor="none",
                    label=f"{display_label} (n={self.tag_value_counts.get(label, 0)})",
                ))
            if omitted_count:
                tag_handles.append(Patch(
                    facecolor="none", edgecolor=self.visual_colors["axis"],
                    label=f"{omitted_count} more values",
                ))
        base_handles = []
        for base in "ACGT":
            base_handles.append(Patch(
                facecolor=self.base_colors[base], edgecolor="none", label=base
            ))
        modification_handles = []
        for label in sorted(self.modification_labels_seen):
            modification_handles.append(Line2D(
                [0], [0], marker="o", linestyle="none",
                markerfacecolor=modification_color(
                    label, self.modification_colors, self.chromosome_palette
                ),
                markeredgecolor="none", label=label,
            ))

        legend_ax = fig.add_axes([
            plot_left, self.legend_bottom_in / fig_height, plot_right - plot_left,
            self.legend_height_in / fig_height,
        ])
        legend_ax.set_xlim(0, 1)
        legend_ax.set_ylim(0, 1)
        legend_ax.set_axis_off()

        alignment_columns = 3 if (pair_geometry_handles or self.molecule_mode) else 1
        alignment_weight = (
            3.35 if pair_geometry_handles else (2.15 if self.molecule_mode else 1.25)
        )
        groups = [
            (
                "Molecules" if self.molecule_mode else "Alignment",
                alignment_handles, alignment_columns, alignment_weight,
            ),
            ("Read events", event_handles, 1, 1.00),
        ]
        if haplotype_handles:
            groups.append(("Haplotype", haplotype_handles, 2, 1.35))
        elif tag_handles:
            tag_columns = 2 if len(tag_handles) > 1 else 1
            tag_weight = 1.35 if len(tag_handles) <= 4 else 1.85
            groups.append((self.tag_label, tag_handles, tag_columns, tag_weight))
        elif insert_size_handles:
            groups.append(("Insert size", insert_size_handles, 1, 1.05))
        if modification_handles:
            groups.append(("Base modifications", modification_handles, 1, 1.10))
        groups.append(("Base identity", base_handles, 2, 1.00))

        positioned_groups = []
        gap = self.styles["legend_compartment_gap"]
        if self.fig_width >= 9:
            total_weight = 0.0
            for group in groups:
                total_weight += group[3]
            usable_width = 1 - gap * (len(groups) - 1)
            cursor = 0.0
            for group in groups:
                width = usable_width * group[3] / total_weight
                positioned_groups.append(group + (cursor, cursor + width, 0.0, 1.0))
                cursor += width + gap
        else:
            column_count = 2 if self.fig_width >= 6 else 1
            row_count = (len(groups) + column_count - 1) // column_count
            cell_width = (1 - gap * (column_count - 1)) / column_count
            cell_height = (1 - gap * (row_count - 1)) / row_count
            for index, group in enumerate(groups):
                row_index = index // column_count
                column_index = index % column_count
                span_columns = (
                    column_count if index == len(groups) - 1 and len(groups) % column_count
                    else 1
                )
                x0 = column_index * (cell_width + gap)
                x1 = x0 + cell_width * span_columns + gap * (span_columns - 1)
                y1 = 1 - row_index * (cell_height + gap)
                y0 = y1 - cell_height
                positioned_groups.append(group + (x0, x1, y0, y1))

        for group in positioned_groups:
            title, handles, columns, weight, x0, x1, y0, y1 = group
            del weight
            legend_ax.add_patch(Rectangle(
                (x0, y0), x1 - x0, y1 - y0,
                transform=legend_ax.transAxes,
                facecolor=self.visual_colors["legend_background"],
                edgecolor=self.visual_colors["legend_edge"],
                linewidth=0.75, zorder=0,
            ))
            legend = legend_ax.legend(
                handles=handles, title=title, loc="center", ncol=columns,
                fontsize=self.styles["legend_font_size"],
                title_fontsize=self.styles["legend_title_size"], frameon=False,
                bbox_to_anchor=((x0 + x1) / 2, (y0 + y1) / 2),
                bbox_transform=legend_ax.transAxes,
                borderaxespad=0, borderpad=0.25,
                columnspacing=1.2, handletextpad=0.45,
                labelspacing=0.35,
            )
            legend.get_title().set_fontweight("bold")
            legend.get_title().set_color(self.visual_colors["primary_text"])
            legend_ax.add_artist(legend)
        return list(legend_ax.artists)

    def separate_legend_from_plots(self, fig, plot_axes, legends: list) -> None:
        """Enforce a physical gap between plot content and the rendered legend."""
        if not legends:
            return
        fig.canvas.draw()
        canvas_renderer = fig.canvas.get_renderer()
        legend_ax = legends[0].axes
        legend_top = legend_ax.get_window_extent(canvas_renderer).y1
        for legend in legends:
            top = legend.get_window_extent(canvas_renderer).y1
            if top > legend_top:
                legend_top = top
        plot_bottom = None
        for ax in plot_axes:
            if ax is legend_ax:
                continue
            bottom = ax.get_tightbbox(canvas_renderer).y0
            if plot_bottom is None or bottom < plot_bottom:
                plot_bottom = bottom
        required_bottom = legend_top + self.legend_plot_gap_in * fig.dpi
        if plot_bottom >= required_bottom:
            return
        shortfall_fraction = (required_bottom - plot_bottom) / fig.bbox.height
        new_bottom = fig.subplotpars.bottom + shortfall_fraction
        maximum_bottom = fig.subplotpars.top - 0.05
        fig.subplots_adjust(bottom=min(new_bottom, maximum_bottom))
        fig.canvas.draw()

    def render_multi(
        self,
        panels: List[dict],
        chrom: str,
        window_start: int,
        window_end: int,
        reference: Optional[ReferenceWindow],
        out_path: str,
        suptitle: str = "",
        genomic_tracks: Optional[List[LoadedAnnotationTrack]] = None,
        contig_length: Optional[int] = None,
        cytobands: Optional[List[Cytoband]] = None,
        assembly_label: Optional[str] = None,
    ) -> None:
        """Stack several BAMs' snapshots in one figure, sharing one genomic
        x-axis - the comparison view for "does aligner A produce
        more gapped alignments than aligner B here". Each panel is a dict with
        keys: label, rows, all_reads_for_coverage, dropped_reads and
        downsampled_reads (optional).
        """
        span = window_end - window_start
        show_ref_track = bool(
            reference and reference.available and
            self.max_reference_span > 0 and span <= self.max_reference_span
        )
        render_base_detail = span <= self.max_mismatch_render_span

        tracks = []
        ratios = []
        if self.show_ideogram and contig_length:
            tracks.append("ideogram")
            ratios.append(self.styles["ideogram_height_in"])
        if show_ref_track:
            tracks.append("reference")
            ratios.append(self.styles["reference_height_in"])
        genomic_tracks = genomic_tracks or []
        for index, annotation in enumerate(genomic_tracks):
            tracks.append(f"annotation_{index}")
            ratios.append(self.annotation_track_height(annotation))

        panel_track_names = []
        coverage_axes = []
        shared_coverage_max = 1
        for i, panel in enumerate(panels):
            n_rows = max(len(panel["rows"]), 1)
            header_name = f"panel_header_{i}"
            cov_name = f"coverage_{i}"
            mod_name = f"modifications_{i}"
            sashimi_name = f"sashimi_{i}"
            aln_name = f"alignments_{i}"
            tracks.append(header_name)
            ratios.append(self.styles["panel_header_height_in"])
            companion_names = []
            for companion_index, annotation in enumerate(panel.get("companion_tracks", [])):
                companion_name = f"companion_{i}_{companion_index}"
                companion_names.append(companion_name)
                tracks.append(companion_name)
                ratios.append(self.annotation_track_height(annotation))
            panel_track_names.append(
                (header_name, cov_name, mod_name, sashimi_name, aln_name, companion_names)
            )
            if self.show_coverage:
                tracks.append(cov_name)
                ratios.append(self.styles["coverage_track_height_in"])
            panel_reads = panel.get("all_reads_for_coverage")
            if panel_reads is None:
                panel_reads = [read for row in panel["rows"] for read in row]
            if (
                self.show_base_modifications
                and has_base_modifications(
                    panel_reads, window_start, window_end, self.modification_codes
                )
            ):
                tracks.append(mod_name)
                ratios.append(self.styles["modification_track_height_in"])
            if self.show_rna_evidence:
                tracks.append(sashimi_name)
                ratios.append(self.styles["sashimi_track_height_in"])
            tracks.append(aln_name)
            ratios.append(max(n_rows * self.row_height_in, self.row_height_in))

        top_margin_in = 0.72
        bottom_margin_in = self.legend_margin_in
        fig_height = sum(ratios) + top_margin_in + bottom_margin_in
        fig, axes = plt.subplots(
            nrows=len(tracks), ncols=1, figsize=(self.fig_width, fig_height), dpi=self.dpi,
            gridspec_kw={"height_ratios": ratios, "hspace": 0.2}, sharex=True,
        )
        ax_by_track = dict(zip(tracks, axes))
        for ax in axes:
            ax.set_xlim(window_start, window_end)
            for spine in ("top", "right", "left"):
                ax.spines[spine].set_visible(False)
            ax.spines["bottom"].set_visible(False)
            ax.tick_params(
                left=False, labelleft=False, bottom=False, top=False,
                labelbottom=False, labeltop=False,
            )

        tick_positions = nice_tick_positions(window_start, window_end)

        title_x, title_ha = self.figure_title_position()
        fig.text(
            title_x, 1 - 0.06 / fig_height,
            f"{chrom}:{window_start + 1:,}-{window_end:,} ({span:,} bp)",
            fontsize=10.5, color=self.visual_colors["primary_text"], fontweight="bold", va="top", ha=title_ha,
        )
        if suptitle:
            fig.text(
                title_x, 1 - 0.34 / fig_height,
                ellipsize(suptitle, max(30, int(self.fig_width * 15))),
                fontsize=8.5, color=self.visual_colors["secondary_text"], va="top", ha=title_ha,
            )

        if "ideogram" in ax_by_track:
            self.draw_ideogram(
                ax_by_track["ideogram"], chrom, window_start, window_end, contig_length,
                cytobands,
            )

        if show_ref_track:
            self.draw_reference_track(
                ax_by_track["reference"], reference, window_start, window_end,
                available_width_in=self.fig_width,
            )

        for index, annotation in enumerate(genomic_tracks):
            self.draw_annotation_track(
                ax_by_track[f"annotation_{index}"], annotation, window_start, window_end
            )

        for track in tracks:
            if track != "ideogram" and not track.startswith("panel_header_"):
                self.draw_background_grid(
                    ax_by_track[track], tick_positions, window_start, window_end
                )
                self.draw_highlights(
                    ax_by_track[track], chrom, window_start, window_end
                )

        for i, panel in enumerate(panels):
            rows = panel["rows"]
            layout = panel.get("layout", "pack")
            n_rows = max(len(rows), 1)
            header_name, cov_name, mod_name, sashimi_name, aln_name, companion_names = panel_track_names[i]

            panel_label = panel.get("label", f"bam{i+1}")
            if panel.get("downsampled_reads"):
                panel_label += f"; {panel['downsampled_reads']} downsampled"
            if panel.get("dropped_reads"):
                panel_label += f"; {panel['dropped_reads']} omitted by --max_rows"
            header_ax = ax_by_track[header_name]
            header_ax.set_ylim(0, 1)
            header_ax.text(
                0.0, 0.45,
                ellipsize(panel_label, max(20, int(self.fig_width * 13))),
                transform=header_ax.transAxes,
                fontsize=9, color=self.visual_colors["primary_text"], fontweight="bold",
                ha="left", va="center", clip_on=True,
            )

            companion_tracks = panel.get("companion_tracks", [])
            for companion_index, annotation in enumerate(companion_tracks):
                self.draw_annotation_track(
                    ax_by_track[companion_names[companion_index]], annotation,
                    window_start, window_end,
                )

            if self.show_coverage:
                cov_ax = ax_by_track[cov_name]
                cov_reads = panel.get("all_reads_for_coverage")
                if not cov_reads:
                    cov_reads = []
                    for row in rows:
                        cov_reads.extend(row)
                panel_coverage_max = self.draw_coverage_track(
                    cov_ax, cov_reads, window_start, window_end
                )
                coverage_axes.append(cov_ax)
                shared_coverage_max = max(shared_coverage_max, panel_coverage_max)

            if mod_name in ax_by_track:
                modification_reads = panel.get("all_reads_for_coverage")
                if modification_reads is None:
                    modification_reads = [read for row in rows for read in row]
                self.draw_modification_track(
                    ax_by_track[mod_name], modification_reads,
                    window_start, window_end,
                )

            if self.show_rna_evidence:
                sashimi_reads = panel.get("all_reads_for_coverage")
                if sashimi_reads is None:
                    sashimi_reads = []
                    for row in rows:
                        sashimi_reads.extend(row)
                self.draw_sashimi_track(
                    ax_by_track[sashimi_name], sashimi_reads,
                    window_start, window_end, chrom=chrom, reference=reference,
                    genomic_tracks=genomic_tracks,
                )

            aln_ax = ax_by_track[aln_name]
            aln_ax.set_ylim(n_rows, 0)
            self.draw_haplotype_lanes(aln_ax, rows)
            for row_idx, row in enumerate(rows):
                y0 = row_idx + self.row_margin
                h = 1 - 2 * self.row_margin
                self.draw_alignment_row(
                    aln_ax, row, y0, h, render_base_detail, layout
                )
            if not rows:
                aln_ax.text(0.5, 0.5, "No alignments in this region", transform=aln_ax.transAxes,
                            ha="center", va="center", fontsize=9,
                            color=self.visual_colors["secondary_text"])
        # Identical y limits make depth differences comparable between sample
        # panels instead of independently stretching every coverage profile.
        if len(coverage_axes) > 1:
            for cov_ax in coverage_axes:
                cov_ax.set_ylim(0, shared_coverage_max * 1.15)
                cov_ax.set_yticks([0, shared_coverage_max])
        for track, ax in ax_by_track.items():
            if track != "ideogram" and not track.startswith("panel_header_"):
                self.draw_center_guide(ax, window_start, window_end)
        bottom_aln_ax = ax_by_track[panel_track_names[-1][4]]
        apply_genomic_axis(
            bottom_aln_ax, tick_positions, window_start, window_end, label_size=9,
            color=self.visual_colors["primary_text"],
        )
        bottom_aln_ax.tick_params(
            bottom=True, labelbottom=True, labelsize=9,
            length=3, colors=self.visual_colors["primary_text"],
        )
        margin_tracks = list(genomic_tracks)
        for panel in panels:
            margin_tracks.extend(panel.get("companion_tracks", []))
        plot_left = left_margin_fraction(self.fig_width, margin_tracks)
        if self.has_split_lanes:
            plot_left = max(plot_left, min(1.15 / self.fig_width, 0.25))
        plot_right = 0.92
        fig.subplots_adjust(left=plot_left, right=plot_right,
                            top=1 - top_margin_in / fig_height,
                            bottom=bottom_margin_in / fig_height)
        self.draw_scale_bar(
            fig, fig_height, plot_left, plot_right, span, assembly_label,
            offset_from_top_in=0.59,
        )
        if self.show_legend:
            legends = self.draw_legends(fig, fig_height, plot_left, plot_right)
            self.separate_legend_from_plots(fig, axes, legends)
        fig.savefig(out_path, dpi=self.dpi)
        plt.close(fig)

    def render_loci(
        self,
        panels: List[dict],
        out_path: str,
        suptitle: str = "",
        assembly_label: Optional[str] = None,
    ) -> None:
        """Render two independently scaled genomic loci as adjacent panels.

        Unlike :meth:`render_multi`, which stacks BAMs over one shared locus,
        each panel here supplies its own chromosome, bounds, and reference.
        This is the IGV-like mate view used to inspect both sides of a
        discordant or split-read event.
        """
        if len(panels) != 2:
            raise ValueError("Mate view requires exactly two locus panels.")

        max_rows = 1
        for panel in panels:
            if len(panel["rows"]) > max_rows:
                max_rows = len(panel["rows"])

        show_ref_track = False
        for panel in panels:
            if (
                panel.get("reference") and panel["reference"].available and
                self.max_reference_span > 0 and
                panel["end"] - panel["start"] <= self.max_reference_span
            ):
                show_ref_track = True
                break

        show_modification_track = self.show_base_modifications and any(
            has_base_modifications(
                panel.get("all_reads_for_coverage")
                if panel.get("all_reads_for_coverage") is not None
                else [read for row in panel["rows"] for read in row],
                panel["start"], panel["end"], self.modification_codes,
            )
            for panel in panels
        )

        tracks = ["panel_header"]
        ratios = [self.styles["panel_header_height_in"]]
        show_ideogram = False
        if self.show_ideogram:
            for panel in panels:
                if panel.get("contig_length"):
                    show_ideogram = True
                    break
        if show_ideogram:
            tracks.append("ideogram")
            ratios.append(self.styles["ideogram_height_in"])
        if show_ref_track:
            tracks.append("reference")
            ratios.append(self.styles["reference_height_in"])
        annotation_count = 0
        for panel in panels:
            track_count = len(panel.get("genomic_tracks", []))
            if track_count > annotation_count:
                annotation_count = track_count
        annotation_row_counts = []
        for index in range(annotation_count):
            shared_rows = 1
            for panel in panels:
                if index < len(panel.get("genomic_tracks", [])):
                    row_count = annotation_row_count(
                        panel["genomic_tracks"][index]
                    )
                    if row_count > shared_rows:
                        shared_rows = row_count
            annotation_row_counts.append(shared_rows)
            tracks.append(f"annotation_{index}")
            shared_height = 0.0
            for panel in panels:
                panel_tracks = panel.get("genomic_tracks", [])
                if index < len(panel_tracks):
                    height = self.annotation_track_height(
                        panel_tracks[index], row_count=shared_rows
                    )
                    shared_height = max(shared_height, height)
            ratios.append(shared_height)
        if self.show_coverage:
            tracks.append("coverage")
            ratios.append(self.styles["coverage_track_height_in"])
        if show_modification_track:
            tracks.append("modifications")
            ratios.append(self.styles["modification_track_height_in"])
        if self.show_rna_evidence:
            tracks.append("sashimi")
            ratios.append(self.styles["sashimi_track_height_in"])
        tracks.append("alignments")
        ratios.append(max(max_rows * self.row_height_in, self.row_height_in))

        top_margin_in = 0.5
        bottom_margin_in = self.legend_margin_in
        fig_height = sum(ratios) + top_margin_in + bottom_margin_in
        fig, axes = plt.subplots(
            nrows=len(tracks), ncols=2, squeeze=False,
            figsize=(self.fig_width, fig_height), dpi=self.dpi,
            gridspec_kw={
                "height_ratios": ratios, "hspace": 0.15,
                "wspace": 0.20 if self.has_split_lanes else 0.12,
            },
        )

        if suptitle:
            title_x, title_ha = self.figure_title_position()
            fig.text(title_x, 0.995, suptitle, fontsize=10.5,
                     color=self.visual_colors["primary_text"],
                     fontweight="bold", va="top", ha=title_ha)

        for panel_idx, panel in enumerate(panels):
            chrom = panel["chrom"]
            start = panel["start"]
            end = panel["end"]
            span = end - start
            rows = panel["rows"]
            layout = panel.get("layout", "pack")
            reference = panel.get("reference")
            self.active_sort_base_position = panel.get(
                "sort_base_position", self.sort_base_position
            )
            self.active_sort_reference_base = panel.get(
                "sort_reference_base", self.sort_reference_base
            )
            render_base_detail = span <= self.max_mismatch_render_span
            axes_by_track = {}
            for i, track in enumerate(tracks):
                axes_by_track[track] = axes[i][panel_idx]
            ticks = nice_tick_positions(start, end, target=4)

            for ax in axes_by_track.values():
                ax.set_xlim(start, end)
                for spine in ("top", "right", "left"):
                    ax.spines[spine].set_visible(False)
                ax.spines["bottom"].set_visible(False)
                ax.tick_params(
                    left=False, labelleft=False, bottom=False, top=False,
                    labelbottom=False, labeltop=False,
                )

            label = panel.get("label", "Primary" if panel_idx == 0 else "Mate")
            dropped = panel.get("dropped_reads", 0)
            dropped_label = f"; {dropped} read(s) omitted" if dropped else ""
            downsampled = panel.get("downsampled_reads", 0)
            downsampled_label = f"; {downsampled} downsampled" if downsampled else ""
            header_ax = axes_by_track["panel_header"]
            header_ax.set_ylim(0, 1)
            header_ax.text(
                0.0, 0.72,
                ellipsize(label, max(16, int(self.fig_width * 6.5))),
                transform=header_ax.transAxes,
                ha="left", va="center", fontsize=8.2,
                color=self.visual_colors["primary_text"],
                fontweight="bold", clip_on=True,
            )
            header_ax.text(
                0.0, 0.16,
                f"{chrom}:{start + 1:,}-{end:,} ({span:,} bp)"
                f"{downsampled_label}{dropped_label}",
                transform=header_ax.transAxes, ha="left", va="center",
                fontsize=7, color=self.visual_colors["secondary_text"], clip_on=True,
            )

            for track in tracks:
                if track not in ("panel_header", "ideogram"):
                    self.draw_background_grid(
                        axes_by_track[track], ticks, start, end
                    )
                    self.draw_highlights(
                        axes_by_track[track], chrom, start, end
                    )

            if show_ref_track:
                ref_ax = axes_by_track["reference"]
                ref_ax.set_ylim(0, 1)
                if (
                    reference and reference.available and
                    self.max_reference_span > 0 and span <= self.max_reference_span
                ):
                    self.draw_reference_track(
                        ref_ax, reference, start, end,
                        available_width_in=self.fig_width / 2,
                    )

            if show_ideogram and panel.get("contig_length"):
                self.draw_ideogram(
                    axes_by_track["ideogram"], chrom, start, end, panel["contig_length"],
                    panel.get("cytobands"),
                )

            panel_annotations = panel.get("genomic_tracks", [])
            for index, annotation in enumerate(panel_annotations):
                self.draw_annotation_track(
                    axes_by_track[f"annotation_{index}"], annotation, start, end,
                    shared_row_count=annotation_row_counts[index],
                )

            if self.show_coverage:
                cov_ax = axes_by_track["coverage"]
                cov_reads = panel.get("all_reads_for_coverage")
                if cov_reads is None:
                    cov_reads = []
                    for row in rows:
                        cov_reads.extend(row)
                self.draw_coverage_track(cov_ax, cov_reads, start, end)

            if show_modification_track:
                modification_reads = panel.get("all_reads_for_coverage")
                if modification_reads is None:
                    modification_reads = [read for row in rows for read in row]
                self.draw_modification_track(
                    axes_by_track["modifications"], modification_reads, start, end
                )

            if self.show_rna_evidence:
                sashimi_reads = panel.get("all_reads_for_coverage")
                if sashimi_reads is None:
                    sashimi_reads = []
                    for row in rows:
                        sashimi_reads.extend(row)
                self.draw_sashimi_track(
                    axes_by_track["sashimi"], sashimi_reads, start, end,
                    chrom=chrom, reference=reference,
                    genomic_tracks=panel_annotations,
                )

            aln_ax = axes_by_track["alignments"]
            aln_ax.set_ylim(max_rows, 0)
            self.draw_haplotype_lanes(aln_ax, rows)
            apply_genomic_axis(
                aln_ax, ticks, start, end, label_size=8,
                color=self.visual_colors["primary_text"],
            )
            aln_ax.tick_params(
                bottom=True, labelbottom=True, labelsize=8,
                length=3, colors=self.visual_colors["primary_text"],
            )
            for row_idx, row in enumerate(rows):
                y0 = row_idx + self.row_margin
                h = 1 - 2 * self.row_margin
                self.draw_alignment_row(
                    aln_ax, row, y0, h, render_base_detail, layout
                )
            if not rows:
                aln_ax.text(0.5, 0.5, "No alignments in this region",
                            transform=aln_ax.transAxes, ha="center", va="center",
                            fontsize=9, color=self.visual_colors["secondary_text"])

            for track, ax in axes_by_track.items():
                if track not in ("panel_header", "ideogram"):
                    self.draw_center_guide(ax, start, end)

        all_genomic_tracks = []
        for panel in panels:
            all_genomic_tracks.extend(panel.get("genomic_tracks", []))
        plot_left = left_margin_fraction(self.fig_width, all_genomic_tracks)
        if self.has_split_lanes:
            plot_left = max(plot_left, min(1.15 / self.fig_width, 0.25))
        plot_right = 0.95
        fig.subplots_adjust(left=plot_left, right=plot_right,
                            top=1 - top_margin_in / fig_height,
                            bottom=bottom_margin_in / fig_height)
        for panel_idx, panel in enumerate(panels):
            panel_bounds = axes[0][panel_idx].get_position()
            self.draw_scale_bar(
                fig, fig_height, panel_bounds.x0, panel_bounds.x1,
                panel["end"] - panel["start"],
                assembly_label if panel_idx == len(panels) - 1 else None,
                offset_from_top_in=0.40,
            )
        if self.show_legend:
            legends = self.draw_legends(fig, fig_height, plot_left, plot_right)
            plot_axes = []
            for row in axes:
                plot_axes.extend(row)
            self.separate_legend_from_plots(fig, plot_axes, legends)
        fig.savefig(out_path, dpi=self.dpi)
        plt.close(fig)

    def render_multi_loci(
        self,
        loci: List[dict],
        out_path: str,
        suptitle: str = "",
        assembly_label: Optional[str] = None,
        link_breakpoints: bool = False,
    ) -> None:
        """Render independently scaled loci as columns and samples as rows.

        Each locus supplies its own chromosome, bounds, reference, annotations,
        and a ``samples`` list.  Sample order must be identical across loci.
        This extends mate view into an explicit multi-locus layout while also
        supporting repeated BAM inputs.
        """
        if len(loci) < 2:
            raise ValueError("Multi-locus view requires at least two loci.")
        sample_count = len(loci[0].get("samples", []))
        if sample_count < 1:
            raise ValueError("Multi-locus view requires at least one sample panel.")
        if any(len(locus.get("samples", [])) != sample_count for locus in loci):
            raise ValueError("Every locus must contain the same number of sample panels.")

        show_ref_track = any(
            locus.get("reference") and locus["reference"].available and
            self.max_reference_span > 0 and
            locus["end"] - locus["start"] <= self.max_reference_span
            for locus in loci
        )
        show_ideogram = self.show_ideogram and any(
            locus.get("contig_length") for locus in loci
        )

        tracks = ["locus_header"]
        ratios = [self.styles["panel_header_height_in"]]
        if show_ideogram:
            tracks.append("ideogram")
            ratios.append(self.styles["ideogram_height_in"])
        if show_ref_track:
            tracks.append("reference")
            ratios.append(self.styles["reference_height_in"])

        annotation_count = max(
            (len(locus.get("genomic_tracks", [])) for locus in loci), default=0
        )
        annotation_row_counts = []
        for annotation_index in range(annotation_count):
            shared_rows = max(
                (
                    annotation_row_count(
                        locus["genomic_tracks"][annotation_index]
                    )
                    for locus in loci
                    if annotation_index < len(locus.get("genomic_tracks", []))
                ),
                default=1,
            )
            annotation_row_counts.append(max(shared_rows, 1))
            tracks.append(f"annotation_{annotation_index}")
            ratios.append(max(
                (
                    self.annotation_track_height(
                        locus["genomic_tracks"][annotation_index],
                        row_count=max(shared_rows, 1),
                    )
                    for locus in loci
                    if annotation_index < len(locus.get("genomic_tracks", []))
                ),
                default=self.styles["annotation_row_height_in"],
            ))

        sample_track_names = []
        for sample_index in range(sample_count):
            header_name = f"sample_header_{sample_index}"
            tracks.append(header_name)
            ratios.append(self.styles["panel_header_height_in"])

            companion_count = max(
                (
                    len(locus["samples"][sample_index].get("companion_tracks", []))
                    for locus in loci
                ),
                default=0,
            )
            companion_names = []
            for companion_index in range(companion_count):
                name = f"companion_{sample_index}_{companion_index}"
                companion_names.append(name)
                tracks.append(name)
                ratios.append(max(
                    (
                        self.annotation_track_height(
                            locus["samples"][sample_index]["companion_tracks"][companion_index]
                        )
                        for locus in loci
                        if companion_index < len(
                            locus["samples"][sample_index].get("companion_tracks", [])
                        )
                    ),
                    default=self.styles["annotation_row_height_in"],
                ))

            coverage_name = f"coverage_{sample_index}"
            modification_name = f"modifications_{sample_index}"
            sashimi_name = f"sashimi_{sample_index}"
            alignment_name = f"alignments_{sample_index}"
            if self.show_coverage:
                tracks.append(coverage_name)
                ratios.append(self.styles["coverage_track_height_in"])
            show_sample_modifications = self.show_base_modifications and any(
                has_base_modifications(
                    locus["samples"][sample_index].get("all_reads_for_coverage")
                    if locus["samples"][sample_index].get("all_reads_for_coverage") is not None
                    else [
                        read
                        for row in locus["samples"][sample_index]["rows"]
                        for read in row
                    ],
                    locus["start"], locus["end"], self.modification_codes,
                )
                for locus in loci
            )
            if show_sample_modifications:
                tracks.append(modification_name)
                ratios.append(self.styles["modification_track_height_in"])
            if self.show_rna_evidence:
                tracks.append(sashimi_name)
                ratios.append(self.styles["sashimi_track_height_in"])
            if self.show_alignments:
                tracks.append(alignment_name)
                maximum_rows = max(
                    len(locus["samples"][sample_index]["rows"])
                    for locus in loci
                )
                ratios.append(max(maximum_rows * self.row_height_in, self.row_height_in))
            sample_track_names.append({
                "header": header_name,
                "companions": companion_names,
                "coverage": coverage_name,
                "modifications": modification_name if show_sample_modifications else None,
                "sashimi": sashimi_name,
                "alignments": alignment_name,
            })

        top_margin_in = 0.58
        bottom_margin_in = self.legend_margin_in
        fig_height = sum(ratios) + top_margin_in + bottom_margin_in
        column_count = len(loci)
        fig, axes = plt.subplots(
            nrows=len(tracks), ncols=column_count, squeeze=False,
            figsize=(self.fig_width, fig_height), dpi=self.dpi,
            gridspec_kw={
                "height_ratios": ratios,
                "hspace": 0.16,
                "wspace": 0.22 if self.has_split_lanes else 0.13,
            },
        )

        if suptitle:
            title_x, title_ha = self.figure_title_position()
            fig.text(
                title_x, 0.995, suptitle, fontsize=10.5,
                color=self.visual_colors["primary_text"], fontweight="bold",
                va="top", ha=title_ha,
            )

        for locus_index, locus in enumerate(loci):
            chrom = locus["chrom"]
            start = locus["start"]
            end = locus["end"]
            span = end - start
            reference = locus.get("reference")
            ticks = nice_tick_positions(start, end, target=max(3, 6 - column_count))
            axes_by_track = {
                track: axes[track_index][locus_index]
                for track_index, track in enumerate(tracks)
            }
            for ax in axes_by_track.values():
                ax.set_xlim(start, end)
                for spine in ("top", "right", "left"):
                    ax.spines[spine].set_visible(False)
                ax.spines["bottom"].set_visible(False)
                ax.tick_params(
                    left=False, labelleft=False, bottom=False, top=False,
                    labelbottom=False, labeltop=False,
                )

            header_ax = axes_by_track["locus_header"]
            header_ax.set_ylim(0, 1)
            header_ax.text(
                0.0, 0.72,
                ellipsize(locus.get("label", f"Locus {locus_index + 1}"),
                          max(12, int(self.fig_width * 13 / column_count))),
                transform=header_ax.transAxes, ha="left", va="center",
                fontsize=8.4, color=self.visual_colors["primary_text"],
                fontweight="bold", clip_on=True,
            )
            header_ax.text(
                0.0, 0.16, f"{chrom}:{start + 1:,}-{end:,} ({span:,} bp)",
                transform=header_ax.transAxes, ha="left", va="center",
                fontsize=7, color=self.visual_colors["secondary_text"], clip_on=True,
            )

            for track, ax in axes_by_track.items():
                if track not in ("locus_header", "ideogram") and not track.startswith("sample_header_"):
                    self.draw_background_grid(ax, ticks, start, end)
                    self.draw_highlights(ax, chrom, start, end)

            if show_ideogram and locus.get("contig_length"):
                self.draw_ideogram(
                    axes_by_track["ideogram"], chrom, start, end,
                    locus["contig_length"], locus.get("cytobands"),
                )
            if show_ref_track:
                reference_ax = axes_by_track["reference"]
                reference_ax.set_ylim(0, 1)
                if (
                    reference and reference.available and self.max_reference_span > 0
                    and span <= self.max_reference_span
                ):
                    self.draw_reference_track(
                        reference_ax, reference, start, end,
                        available_width_in=self.fig_width / column_count,
                    )

            for annotation_index, annotation in enumerate(locus.get("genomic_tracks", [])):
                self.draw_annotation_track(
                    axes_by_track[f"annotation_{annotation_index}"], annotation,
                    start, end, shared_row_count=annotation_row_counts[annotation_index],
                )

            for sample_index, sample in enumerate(locus["samples"]):
                names = sample_track_names[sample_index]
                rows = sample["rows"]
                layout = sample.get("layout", "pack")
                self.active_sort_base_position = sample.get(
                    "sort_base_position", self.sort_base_position
                )
                self.active_sort_reference_base = sample.get(
                    "sort_reference_base", self.sort_reference_base
                )
                render_base_detail = span <= self.max_mismatch_render_span

                sample_label = sample.get("label", f"Sample {sample_index + 1}")
                if sample.get("downsampled_reads"):
                    sample_label += f"; {sample['downsampled_reads']} downsampled"
                if sample.get("dropped_reads"):
                    sample_label += f"; {sample['dropped_reads']} omitted"
                sample_header = axes_by_track[names["header"]]
                sample_header.set_ylim(0, 1)
                sample_header.text(
                    0.0, 0.45,
                    ellipsize(sample_label, max(14, int(self.fig_width * 13 / column_count))),
                    transform=sample_header.transAxes, ha="left", va="center",
                    fontsize=8.2, color=self.visual_colors["primary_text"],
                    fontweight="bold", clip_on=True,
                )

                for companion_index, annotation in enumerate(sample.get("companion_tracks", [])):
                    self.draw_annotation_track(
                        axes_by_track[names["companions"][companion_index]],
                        annotation, start, end,
                    )

                if self.show_coverage:
                    coverage_reads = sample.get("all_reads_for_coverage")
                    if coverage_reads is None:
                        coverage_reads = [read for row in rows for read in row]
                    self.draw_coverage_track(
                        axes_by_track[names["coverage"]], coverage_reads, start, end
                    )
                if names["modifications"]:
                    modification_reads = sample.get("all_reads_for_coverage")
                    if modification_reads is None:
                        modification_reads = [read for row in rows for read in row]
                    self.draw_modification_track(
                        axes_by_track[names["modifications"]],
                        modification_reads, start, end,
                    )
                if self.show_rna_evidence:
                    sashimi_reads = sample.get("all_reads_for_coverage")
                    if sashimi_reads is None:
                        sashimi_reads = [read for row in rows for read in row]
                    self.draw_sashimi_track(
                        axes_by_track[names["sashimi"]], sashimi_reads, start, end,
                        chrom=chrom, reference=reference,
                        genomic_tracks=locus.get("genomic_tracks", []),
                    )
                if self.show_alignments:
                    alignment_ax = axes_by_track[names["alignments"]]
                    maximum_rows = max(
                        len(other_locus["samples"][sample_index]["rows"])
                        for other_locus in loci
                    )
                    alignment_ax.set_ylim(max(maximum_rows, 1), 0)
                    self.draw_haplotype_lanes(alignment_ax, rows)
                    for row_index, row in enumerate(rows):
                        self.draw_alignment_row(
                            alignment_ax, row, row_index + self.row_margin,
                            1 - 2 * self.row_margin, render_base_detail, layout,
                        )
                    if not rows:
                        alignment_ax.text(
                            0.5, 0.5, "No alignments in this region",
                            transform=alignment_ax.transAxes, ha="center", va="center",
                            fontsize=8, color=self.visual_colors["secondary_text"],
                        )

            for track, ax in axes_by_track.items():
                if track not in ("locus_header", "ideogram") and not track.startswith("sample_header_"):
                    self.draw_center_guide(ax, start, end)

            if self.show_alignments:
                axis_track = sample_track_names[-1]["alignments"]
            elif self.show_coverage:
                axis_track = sample_track_names[-1]["coverage"]
            elif annotation_count:
                axis_track = f"annotation_{annotation_count - 1}"
            else:
                axis_track = "reference" if show_ref_track else "locus_header"
            axis_ax = axes_by_track[axis_track]
            apply_genomic_axis(
                axis_ax, ticks, start, end, label_size=8,
                color=self.visual_colors["primary_text"],
            )
            axis_ax.tick_params(
                bottom=True, labelbottom=True, labelsize=8, length=3,
                colors=self.visual_colors["primary_text"],
            )

        all_tracks = []
        for locus in loci:
            all_tracks.extend(locus.get("genomic_tracks", []))
            for sample in locus["samples"]:
                all_tracks.extend(sample.get("companion_tracks", []))
        plot_left = left_margin_fraction(self.fig_width, all_tracks)
        if self.has_split_lanes:
            plot_left = max(plot_left, min(1.15 / self.fig_width, 0.25))
        plot_right = 0.95
        fig.subplots_adjust(
            left=plot_left, right=plot_right,
            top=1 - top_margin_in / fig_height,
            bottom=bottom_margin_in / fig_height,
        )

        for locus_index, locus in enumerate(loci):
            panel_bounds = axes[0][locus_index].get_position()
            self.draw_scale_bar(
                fig, fig_height, panel_bounds.x0, panel_bounds.x1,
                locus["end"] - locus["start"],
                assembly_label if locus_index == len(loci) - 1 else None,
                offset_from_top_in=0.46,
            )

        if link_breakpoints:
            link_color = self.visual_colors["breakpoint_link"]
            link_alpha = self.styles["breakpoint_link_alpha"]
            marker_size = self.styles["breakpoint_link_marker_size"]
            for locus_index in range(len(loci) - 1):
                left_header = axes[0][locus_index]
                right_header = axes[0][locus_index + 1]
                connector = ConnectionPatch(
                    xyA=(0.5, 1.02), coordsA=left_header.transAxes,
                    xyB=(0.5, 1.02), coordsB=right_header.transAxes,
                    axesA=left_header, axesB=right_header,
                    arrowstyle="-", connectionstyle="arc3,rad=-0.08",
                    color=link_color,
                    linewidth=self.styles["breakpoint_link_width"],
                    linestyle=self.styles["breakpoint_link_line_style"],
                    alpha=link_alpha, zorder=30, clip_on=False,
                )
                fig.add_artist(connector)
                left_header.plot(
                    0.5, 1.02, marker="o", markersize=marker_size,
                    color=link_color, alpha=link_alpha,
                    transform=left_header.transAxes, clip_on=False, zorder=31,
                )
                right_header.plot(
                    0.5, 1.02, marker="o", markersize=marker_size,
                    color=link_color, alpha=link_alpha,
                    transform=right_header.transAxes, clip_on=False, zorder=31,
                )

        if self.show_legend and self.show_alignments:
            legends = self.draw_legends(fig, fig_height, plot_left, plot_right)
            self.separate_legend_from_plots(
                fig, [ax for row in axes for ax in row], legends
            )
        fig.savefig(out_path, dpi=self.dpi)
        plt.close(fig)

    def draw_haplotype_lanes(self, ax, rows: List[List[AlignedRead]]) -> None:
        """Shade and label contiguous HP or generic BAM-tag lanes."""
        if not self.has_split_lanes or not rows:
            return
        is_haplotype = self.haplotype_view == "split"
        value_attribute = "haplotype" if is_haplotype else "tag_value"
        lanes = []
        lane_start = 0
        lane_value = getattr(rows[0][0], value_attribute, None) if rows[0] else None
        for row_index, row in enumerate(rows[1:], start=1):
            row_value = getattr(row[0], value_attribute, None) if row else None
            if row_value != lane_value:
                lanes.append((lane_start, row_index, lane_value))
                lane_start = row_index
                lane_value = row_value
        lanes.append((lane_start, len(rows), lane_value))

        for index, (start, end, value) in enumerate(lanes):
            if is_haplotype:
                color = haplotype_color(
                    value, self.haplotype_colors, self.chromosome_palette
                )
            else:
                color = tag_color(value, self.tag_colors, self.chromosome_palette)
            if index % 2 == 0:
                ax.axhspan(
                    start, end, facecolor=color,
                    alpha=self.styles[
                        "haplotype_lane_alpha" if is_haplotype else "tag_lane_alpha"
                    ], zorder=0.2,
                )
            if start:
                ax.axhline(
                    start, color=self.visual_colors["legend_edge"],
                    linewidth=0.7, zorder=1,
                )
            if is_haplotype:
                phase_set_values = set()
                for row in rows[start:end]:
                    for read in row:
                        if getattr(read, "phase_set", None) is not None:
                            phase_set_values.add(str(read.phase_set))
                phase_sets = sorted(phase_set_values)
                label = f"HP {value}" if value is not None else "Untagged"
                if len(phase_sets) == 1:
                    label += f" · PS {phase_sets[0]}"
                elif len(phase_sets) > 1:
                    label += f" · {len(phase_sets)} PS"
            else:
                display_value = str(value) if value is not None else "untagged"
                lane_read_count = sum(len(row) for row in rows[start:end])
                label = f"{self.tag_label}={display_value} (n={lane_read_count})"
            ax.text(
                -0.012, (start + end) / 2, ellipsize(label, 32),
                transform=ax.get_yaxis_transform(), ha="right", va="center",
                fontsize=6.5, color=color, fontweight="bold", clip_on=False,
            )

    def draw_alignment_row(
        self, ax, row: List[AlignedRead], y0: float, h: float,
        render_base_detail: bool, layout: str,
    ) -> None:
        """Draw one row, including visible mate links and one gap annotation."""
        pair_link_width = (
            self.styles["squish_pair_link_width"]
            if self.display_mode == "squish"
            else self.styles["pair_link_width"]
        )
        if self.view_as_pairs:
            pair_members = {}
            for read in row:
                if (
                    read.is_paired and not read.is_secondary
                    and not read.is_supplementary and not read.mate_is_unmapped
                    and read.mate_chrom == read.reference_name
                ):
                    pair_members.setdefault(
                        (read.query_name, read.reference_name), []
                    ).append(read)
            for members in pair_members.values():
                ordered_members = sorted(members, key=attrgetter("ref_start"))
                for index in range(0, len(ordered_members) - 1, 2):
                    left, right = ordered_members[index:index + 2]
                    if left.ref_end < right.ref_start:
                        color, alpha = self.read_style(left)
                        ax.plot(
                            [left.ref_end, right.ref_start],
                            [y0 + h / 2, y0 + h / 2],
                            color=color, alpha=max(alpha, 0.55),
                            linewidth=pair_link_width,
                            zorder=1, solid_capstyle="butt",
                        )

        for read in row:
            self.draw_read(ax, read, y0, h, render_base_detail)

        if (
            self.display_mode != "squish"
            and layout == "expand"
            and self.annotate_gap
        ):
            labels = []
            seen_labels = set()
            for read in row:
                label = read.gap_label()
                if label and label not in seen_labels:
                    seen_labels.add(label)
                    labels.append(label)
            if labels:
                ax.text(
                    1.005, y0 + h / 2, " / ".join(labels),
                    transform=ax.get_yaxis_transform(), fontsize=6.5,
                    va="center", ha="left",
                    color=self.visual_colors["secondary_text"], clip_on=False,
                )

    def draw_read(self, ax, read: AlignedRead, y0: float, h: float, render_base_detail: bool) -> None:
        base_fill, alpha = self.read_style(read)
        squished = self.display_mode == "squish"
        axis_span = abs(ax.get_xlim()[1] - ax.get_xlim()[0])
        pixels_per_base = (
            ax.get_window_extent().width / axis_span if axis_span > 0 else 0
        )
        show_softclip_letters = (
            render_base_detail
            and self.display_mode == "expand"
            and pixels_per_base >= self.styles["softclip_base_letter_min_px"]
        )
        edge_width = (
            self.styles["squish_alignment_edge_width"]
            if squished else self.styles["alignment_edge_width"]
        )
        edge_color = (
            "none" if edge_width == 0
            else self.visual_colors["contrast_edge"]
        )

        if not read.blocks:
            ax.add_patch(Rectangle((read.ref_start, y0), max(read.ref_end - read.ref_start, 1), h,
                                    facecolor=base_fill, alpha=alpha,
                                    edgecolor=edge_color, linewidth=edge_width))
            return

        n_blocks = len(read.blocks)
        for i, blk in enumerate(read.blocks):
            if blk.op in ("M", "=", "X"):
                ax.add_patch(Rectangle((blk.ref_pos, y0), blk.length, h, facecolor=base_fill,
                                        alpha=alpha, edgecolor=edge_color,
                                        linewidth=edge_width, zorder=2))
            elif blk.op in ("D", "N"):
                color = (
                    self.visual_colors["deletion"] if blk.op == "D"
                    else self.visual_colors["reference_skip"]
                )
                line_width = (
                    self.styles[
                        "squish_deletion_line_width"
                        if squished else "deletion_line_width"
                    ]
                ) if blk.op == "D" else self.styles[
                    "squish_split_read_line_width"
                    if squished else "split_read_line_width"
                ]
                ax.plot([blk.ref_pos, blk.ref_pos + blk.length], [y0 + h / 2, y0 + h / 2],
                        color=color, linestyle="-", linewidth=line_width,
                        zorder=3, solid_capstyle="butt")
                if self.show_indel_lengths and not squished and blk.length >= 3:
                    ax.text(blk.ref_pos + blk.length / 2, y0 + h / 2, f"{blk.length}",
                            fontsize=5.5, color=color, ha="center", va="bottom", zorder=4)
            elif blk.op == "I":
                width = self.styles["insertion_marker_width_bp"]
                insertion_color = self.visual_colors["insertion"]
                ax.add_patch(Rectangle((blk.ref_pos - width / 2, y0), width, h, facecolor=insertion_color,
                                        edgecolor="none", linewidth=0, zorder=5))
                if (
                    self.display_mode == "expand"
                    and pixels_per_base * width >= self.styles["insertion_symbol_min_px"]
                ):
                    ax.text(
                        blk.ref_pos, y0 + h / 2, "I",
                        fontsize=self.styles["insertion_symbol_size"],
                        color=self.visual_colors["contrast_edge"],
                        fontweight="bold", ha="center", va="center",
                        clip_on=True, zorder=8,
                    )
                if self.show_indel_lengths and not squished and blk.length >= 3:
                    ax.text(blk.ref_pos, y0, f"+{blk.length}", fontsize=5.5, color=insertion_color,
                            ha="center", va="top", zorder=6)
            elif blk.op == "S":
                is_left = i == 0
                is_right = i == n_blocks - 1
                if not (is_left or is_right):
                    continue
                x0 = blk.ref_pos - blk.length if is_left else blk.ref_pos
                # Soft-clipped bases are real query sequence, just unaligned to the
                # reference - draw them attached to the aligned block, each base
                # colored by its own identity (same convention as mismatches),
                # rather than one flat "clip" color.
                if render_base_detail and read.query_sequence:
                    for offset in range(blk.length):
                        cbase = read.query_sequence[blk.query_pos + offset].upper()
                        cell_x = x0 + offset
                        cell_color = (
                            base_fill if show_softclip_letters
                            else self.base_colors.get(cbase, self.base_colors["N"])
                        )
                        ax.add_patch(Rectangle(
                            (cell_x, y0), 1, h,
                            facecolor=cell_color,
                            alpha=alpha, edgecolor="none", zorder=2,
                        ))
                        if show_softclip_letters:
                            ax.text(
                                cell_x + 0.5, y0 + h / 2, cbase,
                                ha="center", va="center",
                                fontsize=self.styles["softclip_base_letter_size"],
                                color=self.base_colors.get(cbase, self.base_colors["N"]),
                                fontweight="bold",
                                clip_on=True, zorder=8,
                            )
                else:
                    ax.add_patch(Rectangle(
                        (x0, y0), blk.length, h,
                        facecolor=self.visual_colors["softclip"], alpha=alpha,
                        edgecolor=edge_color, linewidth=edge_width,
                        zorder=2,
                    ))
            # 'H' hard clips consume no query bases and are not drawn.

        if render_base_detail:
            for rpos, qbase in read.mismatches:
                ax.add_patch(Rectangle((rpos, y0), 1, h, facecolor=self.base_colors.get(qbase, self.base_colors["N"]),
                                        edgecolor="none", zorder=7))

        if (
            self.show_base_modifications
            and pixels_per_base >= 0.75
        ):
            for modification in read.base_modifications:
                if not modification_matches(modification, self.modification_codes):
                    continue
                probability = modification.probability
                if probability is None or probability < self.min_mod_probability:
                    continue
                label = modification.label
                color = modification_color(
                    label, self.modification_colors, self.chromosome_palette
                )
                self.modification_labels_seen.add(label)
                x = modification.ref_position + 0.5
                ax.scatter(
                    [x], [y0 + h / 2],
                    s=(
                        self.styles["modification_read_marker_size"]
                        * (0.45 if squished else 1.0)
                    ), marker="o",
                    facecolors=color, edgecolors=self.visual_colors["contrast_edge"],
                    linewidths=0.25, alpha=max(0.45, probability), zorder=9,
                    clip_on=True,
                )
                if (
                    not squished
                    and pixels_per_base >= self.styles["modification_letter_min_px"]
                ):
                    letter = str(modification.code)
                    if len(letter) > 2:
                        letter = "•"
                    ax.text(
                        x, y0 + h / 2, letter,
                        ha="center", va="center",
                        fontsize=self.styles["modification_letter_size"],
                        color=self.visual_colors["contrast_edge"],
                        fontweight="bold", clip_on=True, zorder=10,
                    )

        sort_position = self.active_sort_base_position
        if sort_position is not None:
            observed = read.base_at(sort_position)
            reference = self.active_sort_reference_base
            if observed in ("A", "C", "G", "T") and observed != reference:
                ax.add_patch(Rectangle(
                    (sort_position, y0), 1, h,
                    facecolor=self.base_colors[observed], edgecolor="none", zorder=8,
                ))

        if self.molecule_mode and not squished:
            family_size = getattr(read, "molecule_family_size", 1)
            duplicate_reads = getattr(read, "molecule_duplicate_reads", 0)
            is_duplex = getattr(read, "molecule_is_duplex", False)
            if family_size > 1 or duplicate_reads or is_duplex:
                parts = [f"{family_size}×"]
                if duplicate_reads:
                    parts.append(f"dup{duplicate_reads}")
                if is_duplex:
                    parts.append("duplex")
                label = " · ".join(parts)
                read_width_px = max(read.ref_end - read.ref_start, 1) * pixels_per_base
                if read_width_px >= self.styles["molecule_label_min_px"]:
                    x = (read.ref_start + read.ref_end) / 2
                    horizontal_alignment = "center"
                    color = self.visual_colors["contrast_edge"]
                else:
                    x = read.ref_end + max(axis_span * 0.002, 0.5)
                    horizontal_alignment = "left"
                    color = self.visual_colors["secondary_text"]
                ax.text(
                    x, y0 + h / 2, label,
                    ha=horizontal_alignment, va="center",
                    fontsize=self.styles["molecule_label_size"],
                    color=color, fontweight="bold", clip_on=True, zorder=11,
                )
