"""Render many regions from BED/VCF and optionally roll them into one
self-contained, single- or multi-sample HTML report.

BED coordinates are already 0-based half-open, exactly the internal
representation used everywhere else in this codebase (unlike --region, which
is 1-based inclusive and goes through ``cli.parse_region`` first) - so no
coordinate conversion happens here.
"""
from __future__ import annotations

import base64
import html
import os
import re
from dataclasses import dataclass
from datetime import datetime, timezone
from typing import Dict, List, Optional

import pysam

from locus_snap.metrics import RegionSummary

BROWSER_VIEWABLE_FORMATS = {"png", "jpg", "jpeg", "webp", "svg"}
MIME_TYPES = {
    "png": "image/png", "jpg": "image/jpeg", "jpeg": "image/jpeg",
    "webp": "image/webp", "svg": "image/svg+xml",
}
VCF_EXTENSIONS = (".vcf", ".vcf.gz", ".bcf")
_NAME_SANITIZE_RE = re.compile(r"[^A-Za-z0-9._-]+")


@dataclass(frozen=True)
class BatchRegion:
    chrom: str
    start: int  # 0-based, inclusive
    end: int    # 0-based, exclusive
    name: str

    @property
    def display(self) -> str:
        return f"{self.chrom}:{self.start + 1}-{self.end}"


@dataclass
class BatchResult:
    region: BatchRegion
    output_path: Optional[str] = None
    summary: Optional[RegionSummary] = None
    summaries: Optional[List[RegionSummary]] = None
    error: Optional[str] = None

    def sample_summaries(self) -> List[RegionSummary]:
        """Return ordered sample summaries while accepting legacy results."""
        if self.summaries is not None:
            return self.summaries
        return [self.summary] if self.summary is not None else []


def looks_like_vcf(path: str) -> bool:
    return path.lower().endswith(VCF_EXTENSIONS)


def _unique_sanitized_name(raw_name: str, seen_names: Dict[str, int]) -> str:
    """Sanitize a region name for use as a filename/HTML anchor and make it
    unique against every name already produced in this parse."""
    name = _NAME_SANITIZE_RE.sub("_", raw_name.strip()) or "region"
    if name in seen_names:
        seen_names[name] += 1
        name = f"{name}_{seen_names[name]}"
    else:
        seen_names[name] = 1
    return name


def parse_bed_regions(path: str, flank: int = 0) -> List[BatchRegion]:
    """Parse a BED3/BED4 file of regions into 0-based half-open BatchRegions.

    Blank lines and ``#``/``track``/``browser`` header lines are skipped, as
    is standard for BED files. An optional 4th column becomes the region's
    name (and output filename stem); duplicate/missing names are made unique.
    """
    if not os.path.isfile(path):
        raise ValueError(f"Cannot find --batch_regions file: {path}")

    regions: List[BatchRegion] = []
    seen_names: dict = {}
    with open(path, "r", encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, 1):
            line = raw_line.strip()
            if not line:
                continue
            first_token = line.split(None, 1)[0].lower()
            if first_token.startswith("#") or first_token in ("track", "browser"):
                continue
            fields = line.split()
            if len(fields) < 3:
                raise ValueError(
                    f"{path}:{line_number}: expected at least 3 BED columns "
                    f"(chrom, start, end), got: {line!r}"
                )
            chrom = fields[0]
            try:
                start, end = int(fields[1]), int(fields[2])
            except ValueError as exc:
                raise ValueError(
                    f"{path}:{line_number}: BED start/end must be integers"
                ) from exc
            if end <= start:
                raise ValueError(
                    f"{path}:{line_number}: BED end ({end}) must be greater than start ({start})"
                )
            start = max(0, start - flank)
            end = end + flank

            raw_name = fields[3] if len(fields) >= 4 and fields[3] else f"{chrom}_{start}_{end}"
            name = _unique_sanitized_name(raw_name, seen_names)

            regions.append(BatchRegion(chrom=chrom, start=start, end=end, name=name))

    if not regions:
        raise ValueError(f"{path}: no regions found (file is empty or all comments).")
    return regions


def parse_vcf_regions(path: str, flank: int = 0) -> List[BatchRegion]:
    """Turn every ALT-bearing record in a VCF/VCF.gz/BCF into a BatchRegion.

    Uses pysam's resolved ``record.start``/``record.stop`` (which already
    accounts for INFO/END on symbolic and structural records) rather than
    re-deriving span from REF/SVLEN by hand. A multi-allelic line still
    becomes exactly one region - splitting per-ALT would turn one clinical
    variant into several images, which is not what a batch review wants.
    The whole file is read once, top to bottom; no index is required.
    """
    if not os.path.isfile(path):
        raise ValueError(f"Cannot find --batch_regions file: {path}")

    regions: List[BatchRegion] = []
    seen_names: Dict[str, int] = {}
    try:
        with pysam.VariantFile(path) as vcf:
            for record in vcf:
                if not record.alts:
                    continue
                start = max(0, record.start - flank)
                end = record.stop + flank
                if record.id and record.id != ".":
                    raw_name = record.id
                else:
                    alt_text = ",".join(allele[:10] for allele in record.alts)
                    raw_name = f"{record.contig}_{record.pos}_{record.ref[:10]}_{alt_text}"
                name = _unique_sanitized_name(raw_name, seen_names)
                regions.append(BatchRegion(chrom=record.contig, start=start, end=end, name=name))
    except (OSError, ValueError) as exc:
        raise ValueError(f"Cannot read --batch_regions VCF '{path}': {exc}") from exc

    if not regions:
        raise ValueError(f"{path}: no variants with an ALT allele found.")
    return regions


def _summary_cells(summary: Optional[RegionSummary]) -> List[str]:
    if summary is None:
        return ["-", "-", "-", "-"]
    return [
        str(summary.n_reads),
        f"{summary.pct_gapped:.1f}%",
        f"{summary.pct_discordant:.1f}%",
        str(summary.n_softclipped),
    ]


def _report_sample_labels(
    results: List[BatchResult], sample_labels: Optional[List[str]]
) -> List[str]:
    labels = list(sample_labels or [])
    sample_count = max(
        [len(result.sample_summaries()) for result in results] + [len(labels), 1]
    )
    for sample_index in range(len(labels), sample_count):
        inferred = None
        for result in results:
            summaries = result.sample_summaries()
            if sample_index < len(summaries) and summaries[sample_index].label:
                inferred = summaries[sample_index].label
                break
        labels.append(inferred or f"Sample {sample_index + 1}")
    return labels


def _delta(value: float, baseline: float, suffix: str = "") -> str:
    difference = value - baseline
    if abs(difference) < 0.05:
        return "0" + suffix
    return f"{difference:+.1f}{suffix}"


def _sample_metrics_table(
    summaries: List[RegionSummary], sample_labels: List[str]
) -> str:
    if not summaries:
        return ""
    baseline = summaries[0]
    rows = []
    for sample_index, label in enumerate(sample_labels):
        if sample_index >= len(summaries):
            cells = ["-", "-", "-", "-", "-", "-", "-", "-"]
        else:
            summary = summaries[sample_index]
            cells = _summary_cells(summary)
            if sample_index == 0:
                cells.extend(["baseline", "baseline", "baseline", "baseline"])
            else:
                cells.extend([
                    f"{summary.n_reads - baseline.n_reads:+d}",
                    _delta(summary.pct_gapped, baseline.pct_gapped, " pp"),
                    _delta(summary.pct_discordant, baseline.pct_discordant, " pp"),
                    f"{summary.n_softclipped - baseline.n_softclipped:+d}",
                ])
        rows.append(
            "<tr>"
            f"<th scope='row'>{html.escape(label)}</th>"
            + "".join(f"<td>{cell}</td>" for cell in cells)
            + "</tr>"
        )
    return (
        "<div class='table-scroll'><table class='sample-metrics'>"
        "<thead><tr><th>Sample</th><th>Reads</th><th>Gapped</th>"
        "<th>Discordant</th><th>Soft-clipped</th><th>&Delta; reads</th>"
        "<th>&Delta; gapped</th><th>&Delta; discordant</th>"
        "<th>&Delta; soft-clipped</th></tr></thead>"
        f"<tbody>{''.join(rows)}</tbody></table></div>"
    )


def write_html_report(
    results: List[BatchResult],
    report_path: str,
    output_format: str,
    sample_labels: Optional[List[str]] = None,
) -> None:
    """Write one self-contained HTML report with every rendered image embedded
    as a base64 data URI, plus an index table and a per-region summary card.
    """
    selected_format = (output_format or "png").lower().lstrip(".")
    if selected_format not in BROWSER_VIEWABLE_FORMATS:
        choices = ", ".join(sorted(BROWSER_VIEWABLE_FORMATS))
        raise ValueError(
            f"--report needs a browser-viewable --output_format ({choices}); "
            f"got {output_format!r}. PDF/TIFF/SVGZ cannot be embedded in an HTML report."
        )
    mime_type = MIME_TYPES[selected_format]
    labels = _report_sample_labels(results, sample_labels)
    multi_sample = len(labels) > 1

    if multi_sample:
        grouped_headers = "".join(
            f"<th class='sample-group' colspan='4'>{html.escape(label)}</th>"
            for label in labels
        )
        metric_headers = "".join(
            "<th>Reads</th><th>Gapped</th><th>Discordant</th><th>Soft-clipped</th>"
            for _ in labels
        )
        index_header = (
            "<tr><th rowspan='2'>Name</th><th rowspan='2'>Region</th>"
            f"{grouped_headers}</tr><tr>{metric_headers}</tr>"
        )
    else:
        index_header = (
            "<tr><th>Name</th><th>Region</th><th>Reads</th><th>Gapped</th>"
            "<th>Discordant</th><th>Soft-clipped</th></tr>"
        )

    n_ok = 0
    n_failed = 0
    index_rows = []
    cards = []
    for result in results:
        region = result.region
        anchor = html.escape(region.name, quote=True)
        if result.error:
            n_failed += 1
            index_rows.append(
                f"<tr class='failed'><td><a href='#{anchor}'>{html.escape(region.name)}</a></td>"
                f"<td>{html.escape(region.display)}</td>"
                f"<td colspan='{4 * len(labels)}'>failed</td></tr>"
            )
            cards.append(
                f"<section class='card failed' id='{anchor}'>"
                f"<h2>{html.escape(region.name)}</h2>"
                f"<p class='region'>{html.escape(region.display)}</p>"
                f"<p class='error'>{html.escape(result.error)}</p>"
                f"</section>"
            )
            continue

        n_ok += 1
        summaries = result.sample_summaries()
        summary_cells = []
        for sample_index in range(len(labels)):
            summary = summaries[sample_index] if sample_index < len(summaries) else None
            summary_cells.extend(_summary_cells(summary))
        index_rows.append(
            f"<tr><td><a href='#{anchor}'>{html.escape(region.name)}</a></td>"
            f"<td>{html.escape(region.display)}</td>"
            + "".join(f"<td>{cell}</td>" for cell in summary_cells)
            + "</tr>"
        )
        with open(result.output_path, "rb") as image_file:
            encoded = base64.b64encode(image_file.read()).decode("ascii")
        metrics = _sample_metrics_table(summaries, labels) if multi_sample else ""
        if summaries:
            reads, gapped, discordant, softclipped = _summary_cells(summaries[0])
            single_summary = (
                f" &middot; {reads} reads, {gapped} gapped, {discordant} discordant, "
                f"{softclipped} soft-clipped"
            ) if not multi_sample else ""
        else:
            single_summary = ""
        cards.append(
            f"<section class='card' id='{anchor}'>"
            f"<h2>{html.escape(region.name)}</h2>"
            f"<p class='region'>{html.escape(region.display)}{single_summary}</p>"
            f"{metrics}"
            f"<img src='data:{mime_type};base64,{encoded}' alt='{anchor}'>"
            f"</section>"
        )

    generated_at = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    document = f"""<!doctype html>
<html>
<head>
<meta charset="utf-8">
<title>LocusSnap batch report</title>
<style>
  body {{ font-family: -apple-system, Helvetica, Arial, sans-serif; margin: 2rem;
         color: #1a1a1a; background: #fafaf8; }}
  h1 {{ margin-bottom: 0.2rem; }}
  .meta {{ color: #666; margin-top: 0; margin-bottom: 1.5rem; }}
  .table-scroll {{ overflow-x: auto; }}
  table {{ border-collapse: collapse; margin-bottom: 2rem; white-space: nowrap; }}
  th, td {{ padding: 0.35rem 0.75rem; border-bottom: 1px solid #ddd; text-align: left; }}
  thead th {{ background: #f2f3f3; }}
  th.sample-group {{ text-align: center; border-left: 2px solid #fff; }}
  tr.failed td {{ color: #a33; }}
  .card {{ border: 1px solid #ddd; border-radius: 6px; padding: 1rem 1.25rem;
           margin-bottom: 1.5rem; background: #fff; }}
  .card.failed {{ background: #fff5f5; border-color: #e0b4b4; }}
  .card h2 {{ margin-top: 0; margin-bottom: 0.2rem; font-size: 1.05rem; }}
  .card .region {{ color: #555; margin-top: 0; font-size: 0.9rem; }}
  .card .error {{ color: #a33; font-family: monospace; white-space: pre-wrap; }}
  .card img {{ max-width: 100%; height: auto; display: block; margin-top: 0.5rem; }}
  .sample-metrics {{ margin: 0.5rem 0 1rem; font-size: 0.84rem; }}
  .sample-metrics th[scope='row'] {{ background: transparent; }}
</style>
</head>
<body>
<h1>LocusSnap batch report</h1>
<p class="meta">Generated {generated_at} &middot; {n_ok} rendered, {n_failed} failed &middot; {len(labels)} sample{'s' if len(labels) != 1 else ''}</p>
<div class="table-scroll"><table>
<thead>{index_header}</thead>
<tbody>
{"".join(index_rows)}
</tbody>
</table></div>
{"".join(cards)}
</body>
</html>
"""
    parent = os.path.dirname(report_path)
    if parent:
        os.makedirs(parent, exist_ok=True)
    with open(report_path, "w", encoding="utf-8") as handle:
        handle.write(document)
