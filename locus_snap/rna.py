"""RNA-seq splice-junction and chimeric-fusion evidence."""
from __future__ import annotations

from dataclasses import dataclass, field
import csv
import re
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

from locus_snap.annotations import LoadedAnnotationTrack
from locus_snap.read_model import AlignedRead, CigarBlock, SAEntry
from locus_snap.reference import ReferenceWindow


RNA_STRANDNESS_MODES = ("alignment", "forward", "reverse", "unstranded")
CANONICAL_SPLICE_MOTIFS = {"GT-AG", "GC-AG", "AT-AC"}
_CIGAR_PATTERN = re.compile(r"(\d+)([MIDNSHP=X])")


@dataclass(frozen=True)
class JunctionEvidence:
    donor: int
    acceptor: int
    strand: str
    read_names: Tuple[str, ...]
    minimum_left_anchor: int
    minimum_right_anchor: int

    @property
    def count(self) -> int:
        return len(self.read_names)


@dataclass(frozen=True)
class FusionEvidence:
    local_breakpoint: int
    partner_chrom: str
    partner_breakpoint: int
    local_strand: str
    partner_strand: str
    split_read_names: Tuple[str, ...] = ()
    spanning_pair_names: Tuple[str, ...] = ()

    @property
    def split_reads(self) -> int:
        return len(self.split_read_names)

    @property
    def spanning_pairs(self) -> int:
        return len(self.spanning_pair_names)

    @property
    def support(self) -> int:
        return len(set(self.split_read_names) | set(self.spanning_pair_names))


@dataclass
class _JunctionAccumulator:
    names: Set[str] = field(default_factory=set)
    left_anchors: List[int] = field(default_factory=list)
    right_anchors: List[int] = field(default_factory=list)


@dataclass
class _FusionAccumulator:
    local_positions: List[int] = field(default_factory=list)
    partner_positions: List[int] = field(default_factory=list)
    local_strands: List[str] = field(default_factory=list)
    partner_strands: List[str] = field(default_factory=list)
    split_names: Set[str] = field(default_factory=set)
    pair_names: Set[str] = field(default_factory=set)


def _opposite_strand(strand: str) -> str:
    return "-" if strand == "+" else "+" if strand == "-" else "."


def transcript_strand(read: AlignedRead, strandness: str = "alignment") -> str:
    """Map an alignment strand to the inferred transcript strand.

    ``forward`` and ``reverse`` follow common paired-end RNA library rules:
    read 1 defines orientation and read 2 is normalized to the opposite mate.
    ``alignment`` preserves the historical LocusSnap behaviour.
    """
    if strandness not in RNA_STRANDNESS_MODES:
        choices = ", ".join(RNA_STRANDNESS_MODES)
        raise ValueError(f"RNA strandness must be one of: {choices}.")
    if strandness == "unstranded":
        return "."
    strand = getattr(read, "strand", ".")
    if strandness == "alignment":
        return strand
    if getattr(read, "is_read2", False):
        strand = _opposite_strand(strand)
    if strandness == "reverse":
        strand = _opposite_strand(strand)
    return strand


def _junction_blocks(read: AlignedRead):
    blocks = list(getattr(read, "blocks", []) or [])
    for index, block in enumerate(blocks):
        if block.op != "N" or block.length <= 0:
            continue
        left_anchor = 0
        right_anchor = 0
        for candidate in reversed(blocks[:index]):
            if candidate.op in ("M", "=", "X"):
                left_anchor = candidate.length
                break
            if candidate.op in ("N", "D"):
                break
        for candidate in blocks[index + 1:]:
            if candidate.op in ("M", "=", "X"):
                right_anchor = candidate.length
                break
            if candidate.op in ("N", "D"):
                break
        yield block.ref_pos, block.ref_pos + block.length, left_anchor, right_anchor


def collect_splice_junctions(
    reads: Sequence[AlignedRead], strand_mode: str = "combined",
    strandness: str = "alignment", minimum_anchor: int = 0,
) -> Dict[Tuple[int, int, str], JunctionEvidence]:
    """Collect unique-read support for CIGAR-N splice junctions."""
    if strand_mode not in ("combined", "split"):
        raise ValueError("Sashimi strand mode must be combined or split.")
    if minimum_anchor < 0:
        raise ValueError("Minimum junction anchor cannot be negative.")
    accumulators: Dict[Tuple[int, int, str], _JunctionAccumulator] = {}
    for read_index, read in enumerate(reads):
        strand = transcript_strand(read, strandness) if strand_mode == "split" else "."
        name = str(getattr(read, "query_name", f"read-{read_index}"))
        observed = list(_junction_blocks(read))
        if not observed:
            # Compatibility for lightweight callers that expose only the
            # historical deletion tuple representation.
            observed = [
                (position, position + length, minimum_anchor, minimum_anchor)
                for position, length, is_skip in getattr(read, "deletions", [])
                if is_skip and length > 0
            ]
        for donor, acceptor, left_anchor, right_anchor in observed:
            if left_anchor < minimum_anchor or right_anchor < minimum_anchor:
                continue
            key = (donor, acceptor, strand)
            accumulator = accumulators.setdefault(key, _JunctionAccumulator())
            accumulator.names.add(name)
            accumulator.left_anchors.append(left_anchor)
            accumulator.right_anchors.append(right_anchor)
    return {
        key: JunctionEvidence(
            donor=key[0], acceptor=key[1], strand=key[2],
            read_names=tuple(sorted(accumulator.names)),
            minimum_left_anchor=min(accumulator.left_anchors),
            minimum_right_anchor=min(accumulator.right_anchors),
        )
        for key, accumulator in accumulators.items()
    }


def annotated_junctions(
    tracks: Optional[Sequence[LoadedAnnotationTrack]],
) -> Set[Tuple[str, int, int, str]]:
    """Extract exon-boundary junctions from visible BED12/GFF/GTF models."""
    known: Set[Tuple[str, int, int, str]] = set()
    for track in tracks or []:
        if track.kind not in ("bed", "gff", "gff3", "gtf"):
            continue
        chrom = track.chrom
        for item in track.items:
            intervals = sorted(set(list(item.blocks) + list(item.utrs)))
            merged: List[List[int]] = []
            for start, end in intervals:
                if end <= start:
                    continue
                if not merged or start > merged[-1][1]:
                    merged.append([start, end])
                else:
                    merged[-1][1] = max(merged[-1][1], end)
            for left, right in zip(merged, merged[1:]):
                if left[1] < right[0]:
                    known.add((chrom, left[1], right[0], item.strand))
    return known


def is_annotated_junction(
    chrom: str, junction: JunctionEvidence,
    known: Set[Tuple[str, int, int, str]],
) -> bool:
    for strand in (junction.strand, "."):
        if (chrom, junction.donor, junction.acceptor, strand) in known:
            return True
    if junction.strand == ".":
        return any(
            (chrom, junction.donor, junction.acceptor, strand) in known
            for strand in ("+", "-")
        )
    return False


def splice_motif(
    junction: JunctionEvidence, reference: Optional[ReferenceWindow],
) -> Optional[str]:
    """Return the strand-oriented donor-acceptor dinucleotide motif."""
    if reference is None or not reference.available:
        return None
    donor = "".join(
        reference.base_at(junction.donor + offset) or "N" for offset in range(2)
    )
    acceptor = "".join(
        reference.base_at(junction.acceptor - 2 + offset) or "N" for offset in range(2)
    )
    if "N" in donor + acceptor:
        return None
    plus = f"{donor}-{acceptor}"
    if junction.strand != "-":
        if junction.strand == "." and plus not in CANONICAL_SPLICE_MOTIFS:
            reverse = _reverse_complement(acceptor) + "-" + _reverse_complement(donor)
            return reverse if reverse in CANONICAL_SPLICE_MOTIFS else plus
        return plus
    return _reverse_complement(acceptor) + "-" + _reverse_complement(donor)


def _reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def _terminal_clips(cigar: str) -> Tuple[int, int]:
    operations = _CIGAR_PATTERN.findall(cigar or "")
    if not operations:
        return 0, 0
    left = int(operations[0][0]) if operations[0][1] in ("S", "H") else 0
    right = int(operations[-1][0]) if operations[-1][1] in ("S", "H") else 0
    return left, right


def _local_breakpoint(read: AlignedRead) -> int:
    left = int(getattr(read, "soft_clip_left", 0)) + int(
        getattr(read, "hard_clip_left", 0)
    )
    right = int(getattr(read, "soft_clip_right", 0)) + int(
        getattr(read, "hard_clip_right", 0)
    )
    if left > right:
        return int(read.ref_start)
    if right > left:
        return int(read.ref_end)
    return int(read.ref_start if getattr(read, "is_reverse", False) else read.ref_end)


def _sa_breakpoint(entry: SAEntry) -> int:
    left, right = _terminal_clips(entry.cigar)
    if left > right:
        return entry.start
    if right > left:
        return entry.end
    return entry.start if entry.strand == "+" else entry.end


def _mate_breakpoint(read: AlignedRead, mate_start: int) -> int:
    if getattr(read, "mate_is_reverse", False):
        return mate_start
    aligned_span = max(int(getattr(read, "ref_end", 0)) - int(getattr(read, "ref_start", 0)), 1)
    return mate_start + aligned_span


def _median(values: Sequence[int]) -> int:
    ordered = sorted(values)
    return ordered[len(ordered) // 2]


def _majority_strand(values: Iterable[str]) -> str:
    counts = {"+": 0, "-": 0, ".": 0}
    for value in values:
        counts[value if value in counts else "."] += 1
    return max(counts, key=lambda value: (counts[value], value))


def collect_fusion_evidence(
    reads: Sequence[AlignedRead], chrom: str, breakpoint_tolerance: int = 10,
    minimum_distance: int = 100_000, minimum_mapq: int = 20,
) -> List[FusionEvidence]:
    """Cluster SA split reads and spanning pairs into candidate fusions."""
    if breakpoint_tolerance < 0:
        raise ValueError("Fusion breakpoint tolerance cannot be negative.")
    if minimum_distance < 0:
        raise ValueError("Fusion minimum distance cannot be negative.")
    if minimum_mapq < 0:
        raise ValueError("Fusion minimum MAPQ cannot be negative.")
    clusters: List[Tuple[str, _FusionAccumulator]] = []

    def add_event(
        local: int, partner_chrom: str, partner: int,
        local_strand: str, partner_strand: str, name: str, evidence_type: str,
    ) -> None:
        selected = None
        for candidate_chrom, cluster in clusters:
            if candidate_chrom != partner_chrom:
                continue
            if (
                abs(local - _median(cluster.local_positions)) <= breakpoint_tolerance
                and abs(partner - _median(cluster.partner_positions)) <= breakpoint_tolerance
            ):
                selected = cluster
                break
        if selected is None:
            selected = _FusionAccumulator()
            clusters.append((partner_chrom, selected))
        selected.local_positions.append(local)
        selected.partner_positions.append(partner)
        selected.local_strands.append(local_strand)
        selected.partner_strands.append(partner_strand)
        if evidence_type == "split":
            selected.split_names.add(name)
        else:
            selected.pair_names.add(name)

    for read_index, read in enumerate(reads):
        if int(getattr(read, "mapq", 0)) < minimum_mapq:
            continue
        name = str(getattr(read, "query_name", f"read-{read_index}"))
        local = _local_breakpoint(read)
        local_strand = getattr(read, "strand", ".")
        for entry in getattr(read, "sa_entries", []) or []:
            if entry.mapq < minimum_mapq:
                continue
            distant = entry.rname != chrom or abs(_sa_breakpoint(entry) - local) >= minimum_distance
            if not distant:
                continue
            add_event(
                local, entry.rname, _sa_breakpoint(entry), local_strand,
                entry.strand, name, "split",
            )

        mate_chrom = getattr(read, "mate_chrom", None)
        mate_start = getattr(read, "mate_start", None)
        if (
            mate_chrom is not None and mate_start is not None
            and not getattr(read, "mate_is_unmapped", False)
        ):
            partner_breakpoint = _mate_breakpoint(read, int(mate_start))
            distant = mate_chrom != chrom or abs(partner_breakpoint - local) >= minimum_distance
            if distant:
                add_event(
                    local, str(mate_chrom), partner_breakpoint, local_strand, ".",
                    name, "pair",
                )

    evidence = []
    for partner_chrom, cluster in clusters:
        evidence.append(FusionEvidence(
            local_breakpoint=_median(cluster.local_positions),
            partner_chrom=partner_chrom,
            partner_breakpoint=_median(cluster.partner_positions),
            local_strand=_majority_strand(cluster.local_strands),
            partner_strand=_majority_strand(cluster.partner_strands),
            split_read_names=tuple(sorted(cluster.split_names)),
            spanning_pair_names=tuple(sorted(cluster.pair_names)),
        ))
    return sorted(
        evidence,
        key=lambda item: (-item.support, item.local_breakpoint, item.partner_chrom),
    )


def gene_label_at(
    tracks: Optional[Sequence[LoadedAnnotationTrack]], chrom: str, position: int,
) -> Optional[str]:
    """Return the narrowest visible gene/transcript label covering a position."""
    candidates = []
    for track in tracks or []:
        if track.chrom != chrom or track.kind not in ("bed", "gff", "gff3", "gtf"):
            continue
        for item in track.items:
            if item.start <= position < item.end:
                label = item.group_label or item.name or item.transcript_label
                if label:
                    candidates.append((item.end - item.start, label))
    return min(candidates)[1] if candidates else None


RNA_EVIDENCE_FIELDS = (
    "event_type", "chrom", "start", "end", "strand", "status", "motif",
    "support", "minimum_left_anchor", "minimum_right_anchor",
    "partner_chrom", "partner_position", "split_reads", "spanning_pairs",
    "read_names",
)


def rna_evidence_rows(
    reads: Sequence[AlignedRead], chrom: str,
    reference: Optional[ReferenceWindow] = None,
    genomic_tracks: Optional[Sequence[LoadedAnnotationTrack]] = None,
    strand_mode: str = "combined", strandness: str = "alignment",
    minimum_junction_reads: int = 1, minimum_anchor: int = 0,
    include_fusions: bool = True, minimum_fusion_reads: int = 2,
    fusion_breakpoint_tolerance: int = 10,
    fusion_minimum_distance: int = 100_000, minimum_fusion_mapq: int = 20,
) -> List[dict]:
    """Return stable, tabular junction and fusion summaries."""
    known = annotated_junctions(genomic_tracks)
    annotation_available = any(
        track.kind in ("bed", "gff", "gff3", "gtf") and track.items
        for track in genomic_tracks or []
    )
    rows = []
    junctions = collect_splice_junctions(
        reads, strand_mode=strand_mode, strandness=strandness,
        minimum_anchor=minimum_anchor,
    )
    for evidence in sorted(
        junctions.values(), key=lambda item: (item.donor, item.acceptor, item.strand)
    ):
        if evidence.count < minimum_junction_reads:
            continue
        annotated = is_annotated_junction(chrom, evidence, known)
        status = "annotated" if annotated else "novel" if annotation_available else "unclassified"
        rows.append({
            "event_type": "junction", "chrom": chrom,
            "start": evidence.donor + 1, "end": evidence.acceptor,
            "strand": evidence.strand, "status": status,
            "motif": splice_motif(evidence, reference) or "",
            "support": evidence.count,
            "minimum_left_anchor": evidence.minimum_left_anchor,
            "minimum_right_anchor": evidence.minimum_right_anchor,
            "partner_chrom": "", "partner_position": "",
            "split_reads": "", "spanning_pairs": "",
            "read_names": ",".join(evidence.read_names),
        })
    if include_fusions:
        for evidence in collect_fusion_evidence(
            reads, chrom, breakpoint_tolerance=fusion_breakpoint_tolerance,
            minimum_distance=fusion_minimum_distance,
            minimum_mapq=minimum_fusion_mapq,
        ):
            if evidence.support < minimum_fusion_reads:
                continue
            rows.append({
                "event_type": "fusion", "chrom": chrom,
                "start": evidence.local_breakpoint + 1,
                "end": evidence.local_breakpoint + 1,
                "strand": evidence.local_strand, "status": "candidate",
                "motif": "", "support": evidence.support,
                "minimum_left_anchor": "", "minimum_right_anchor": "",
                "partner_chrom": evidence.partner_chrom,
                "partner_position": evidence.partner_breakpoint + 1,
                "split_reads": evidence.split_reads,
                "spanning_pairs": evidence.spanning_pairs,
                "read_names": ",".join(sorted(
                    set(evidence.split_read_names) | set(evidence.spanning_pair_names)
                )),
            })
    return rows


def write_rna_evidence_tsv(
    path: str, reads: Sequence[AlignedRead], chrom: str, **kwargs,
) -> None:
    """Write one row per supported splice junction or fusion candidate."""
    rows = rna_evidence_rows(reads, chrom, **kwargs)
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=RNA_EVIDENCE_FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
