"""Build molecule-family consensus alignments from standard SAM tags."""
from __future__ import annotations

from collections import Counter, defaultdict
from copy import copy
from dataclasses import dataclass
from statistics import median
from typing import Dict, List, Optional, Sequence, Tuple

from locus_snap.read_model import AlignedRead, CigarBlock
from locus_snap.reference import ReferenceWindow


MOLECULE_TAGS = ("MI", "RX", "UB")


@dataclass(frozen=True)
class MoleculeBuildResult:
    reads: List[AlignedRead]
    tag: str
    input_read_count: int
    tagged_read_count: int
    family_count: int
    filtered_family_count: int


def resolve_molecule_tag(reads: Sequence[AlignedRead], requested: str = "auto") -> str:
    """Resolve an explicit tag or choose the best-covered standard tag."""
    selected = requested.upper()
    if selected != "AUTO":
        if selected not in MOLECULE_TAGS:
            choices = ", ".join(MOLECULE_TAGS)
            raise ValueError(f"Molecule tag must be auto or one of: {choices}.")
        if reads and not any(selected in read.molecule_tags for read in reads):
            raise ValueError(f"No fetched reads carry molecule tag {selected}.")
        return selected
    if not reads:
        return "MI"
    counts = {
        tag: sum(tag in read.molecule_tags for read in reads)
        for tag in MOLECULE_TAGS
    }
    best = max(MOLECULE_TAGS, key=lambda tag: (counts[tag], -MOLECULE_TAGS.index(tag)))
    if counts[best] == 0:
        raise ValueError(
            "Molecule mode could not find MI, RX, or UB tags in the fetched reads."
        )
    return best


def _read_side(read: AlignedRead) -> int:
    if read.is_read1:
        return 1
    if read.is_read2:
        return 2
    return 0


def _family_groups(
    reads: Sequence[AlignedRead], tag: str, position_tolerance: int,
) -> List[List[AlignedRead]]:
    grouped: Dict[Tuple[str, str, str, int], List[List[AlignedRead]]] = defaultdict(list)
    singleton_families = []
    ordered = sorted(
        reads,
        key=lambda read: (
            read.reference_name, read.ref_start, read.ref_end, read.query_name,
        ),
    )
    for read in ordered:
        value = read.molecule_tags.get(tag)
        if value is None:
            singleton_families.append([read])
            continue
        cell = read.molecule_tags.get("CB", "")
        key = (value, cell, read.reference_name, _read_side(read))
        candidates = grouped[key]
        selected_family = None
        for family in candidates:
            anchor = family[0]
            if (
                abs(read.ref_start - anchor.ref_start) <= position_tolerance
                and abs(read.ref_end - anchor.ref_end) <= position_tolerance
            ):
                selected_family = family
                break
        if selected_family is None:
            selected_family = []
            candidates.append(selected_family)
        selected_family.append(read)

    families = list(singleton_families)
    for candidates in grouped.values():
        families.extend(candidates)
    return sorted(
        families,
        key=lambda family: (
            min(read.ref_start for read in family),
            min(read.query_name for read in family),
        ),
    )


def _base_quality(read: AlignedRead, position: int) -> int:
    for block in read.blocks:
        if block.op not in ("M", "=", "X"):
            continue
        if block.ref_pos <= position < block.ref_pos + block.length:
            query_index = block.query_pos + position - block.ref_pos
            if query_index < len(read.query_qualities):
                return int(read.query_qualities[query_index])
            return 0
    return 0


def _majority(values: Sequence[Optional[str]]) -> Optional[str]:
    counts = Counter(value for value in values if value is not None)
    if not counts:
        return None
    return max(counts, key=lambda value: (counts[value], str(value)))


def _supported_events(
    families: Sequence[AlignedRead], attribute: str, minimum_fraction: float,
):
    counts = Counter()
    for read in families:
        counts.update(set(getattr(read, attribute)))
    threshold = minimum_fraction * len(families)
    return sorted(event for event, count in counts.items() if count >= threshold)


def build_consensus_read(
    family: Sequence[AlignedRead], tag: str,
    minimum_consensus_fraction: float = 0.60,
    reference: Optional[ReferenceWindow] = None,
) -> AlignedRead:
    """Collapse one positional UMI family into a drawable consensus alignment."""
    representative = max(
        family,
        key=lambda read: (
            not read.is_duplicate, read.mapq, -read.mismatch_count, read.query_name,
        ),
    )
    consensus = copy(representative)
    family_size = len(family)
    start = min(read.ref_start for read in family)
    end = max(read.ref_end for read in family)
    calls: Dict[int, str] = {}
    qualities: Dict[int, int] = {}
    fractions = []
    for position in range(start, end):
        observations = []
        quality_sums = Counter()
        for read in family:
            base = read.base_at(position)
            if base not in ("A", "C", "G", "T", "-", "~"):
                continue
            observations.append(base)
            if base in "ACGT":
                quality_sums[base] += _base_quality(read, position)
        if not observations:
            continue
        counts = Counter(observations)
        winning_base = max(
            counts,
            key=lambda base: (counts[base], quality_sums[base], base),
        )
        # Missing coverage is not a vote: a one-read overhang must not extend
        # an otherwise large molecule family as a high-confidence consensus.
        fraction = counts[winning_base] / family_size
        if fraction < minimum_consensus_fraction:
            continue
        calls[position] = winning_base
        fractions.append(fraction)
        if winning_base in "ACGT":
            mean_quality = quality_sums[winning_base] / counts[winning_base]
            qualities[position] = min(60, int(round(mean_quality * fraction)))

    blocks = []
    query_sequence = []
    query_qualities = []
    positions = sorted(calls)
    index = 0
    while index < len(positions):
        run_start = positions[index]
        base = calls[run_start]
        op = "M" if base in "ACGT" else ("D" if base == "-" else "N")
        run_positions = [run_start]
        index += 1
        while index < len(positions):
            position = positions[index]
            next_base = calls[position]
            next_op = "M" if next_base in "ACGT" else ("D" if next_base == "-" else "N")
            if position != run_positions[-1] + 1 or next_op != op:
                break
            run_positions.append(position)
            index += 1
        query_position = len(query_sequence)
        if op == "M":
            for position in run_positions:
                query_sequence.append(calls[position])
                query_qualities.append(qualities.get(position, 0))
        blocks.append(CigarBlock(op, run_start, query_position, len(run_positions)))

    insertions = _supported_events(family, "insertions", minimum_consensus_fraction)
    for position, length in insertions:
        blocks.append(CigarBlock("I", position, len(query_sequence), length))
    blocks.sort(key=lambda block: (block.ref_pos, block.op != "I"))

    deletions = []
    for block in blocks:
        if block.op in ("D", "N"):
            deletions.append((block.ref_pos, block.length, block.op == "N"))

    consensus.ref_start = start
    consensus.ref_end = end
    consensus.blocks = blocks
    consensus.query_sequence = "".join(query_sequence)
    consensus.query_qualities = query_qualities
    consensus.consensus_bases = calls
    consensus.consensus_qualities = qualities
    consensus.insertions = insertions
    consensus.deletions = deletions
    consensus.soft_clip_left = 0
    consensus.soft_clip_right = 0
    consensus.soft_clip_total = 0
    consensus.hard_clip_left = 0
    consensus.hard_clip_right = 0
    consensus.base_modifications = []

    mismatches = []
    mismatch_details = []
    if reference is not None and reference.available:
        for position, base in calls.items():
            reference_base = reference.base_at(position)
            if (
                base in "ACGT"
                and reference_base is not None
                and reference_base in "ACGT"
                and base != reference_base
            ):
                mismatches.append((position, base))
                mismatch_details.append((position, base, qualities.get(position, 0)))
    consensus.mismatches = mismatches
    consensus.mismatch_details = mismatch_details
    consensus.mismatch_count = len(mismatches)

    deletion_length = sum(length for _, length, is_skip in deletions if not is_skip)
    consensus.cigar_gap_len = deletion_length + sum(length for _, length in insertions)
    consensus.sa_gap_len = max(read.sa_gap_len for read in family)
    consensus.gap_length = max(consensus.cigar_gap_len, consensus.sa_gap_len)
    consensus.sa_count = max(read.sa_count for read in family)
    consensus.has_cross_chrom_sa = any(read.has_cross_chrom_sa for read in family)
    consensus.mapq = int(round(sum(read.mapq for read in family) / family_size))
    consensus.insert_size = int(round(median(read.insert_size for read in family)))
    consensus.haplotype = _majority([read.haplotype for read in family])
    consensus.phase_set = _majority([read.phase_set for read in family])
    consensus.tag_value = _majority([read.tag_value for read in family])
    consensus.is_duplicate = False
    consensus.is_secondary = all(read.is_secondary for read in family)
    consensus.is_supplementary = all(read.is_supplementary for read in family)
    consensus.molecule_tag = tag
    molecule_value = representative.molecule_tags.get(tag, "untagged")
    cell = representative.molecule_tags.get("CB")
    molecule_prefix = f"{cell}:" if cell else ""
    consensus.molecule_id = (
        f"{tag}:{molecule_prefix}{molecule_value}:"
        f"{representative.reference_name}:{start + 1}-{end}"
    )
    consensus.query_name = consensus.molecule_id
    consensus.molecule_family_size = family_size
    consensus.molecule_duplicate_reads = sum(read.is_duplicate for read in family)
    consensus.molecule_member_names = tuple(sorted({read.query_name for read in family}))
    consensus.molecule_is_duplex = len({read.strand for read in family}) > 1
    consensus.molecule_consensus_fraction = (
        sum(fractions) / len(fractions) if fractions else 0.0
    )
    consensus.strand = representative.strand
    consensus.is_reverse = representative.is_reverse
    consensus.pair_category = _majority([read.pair_category for read in family]) or "normal"
    return consensus


def build_molecule_consensus_reads(
    reads: Sequence[AlignedRead], requested_tag: str = "auto",
    minimum_family_size: int = 1, position_tolerance: int = 2,
    minimum_consensus_fraction: float = 0.60,
    reference: Optional[ReferenceWindow] = None,
) -> MoleculeBuildResult:
    """Group positional UMI families and return one consensus unit per family."""
    if minimum_family_size < 1:
        raise ValueError("Minimum molecule family size must be at least one.")
    if position_tolerance < 0:
        raise ValueError("Molecule position tolerance cannot be negative.")
    if not 0.5 <= minimum_consensus_fraction <= 1:
        raise ValueError("Molecule consensus fraction must be between 0.5 and 1.")
    tag = resolve_molecule_tag(reads, requested_tag)
    families = _family_groups(reads, tag, position_tolerance)
    retained = [family for family in families if len(family) >= minimum_family_size]
    consensus_reads = [
        build_consensus_read(
            family, tag, minimum_consensus_fraction=minimum_consensus_fraction,
            reference=reference,
        )
        for family in retained
    ]
    return MoleculeBuildResult(
        reads=consensus_reads,
        tag=tag,
        input_read_count=len(reads),
        tagged_read_count=sum(tag in read.molecule_tags for read in reads),
        family_count=len(families),
        filtered_family_count=len(families) - len(retained),
    )
