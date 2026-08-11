#!/usr/bin/env python3
"""Regenerate the deterministic BAM and tabix files used by demo figures."""
from dataclasses import dataclass
from array import array
from math import exp, log2
from pathlib import Path
from random import Random

import pysam


PROJECT_DIR = Path(__file__).resolve().parent
OUT_DIR = PROJECT_DIR / "out"
DEMO_DATA_DIR = OUT_DIR / "demo_data"
ALIGNMENTS_DIR = DEMO_DATA_DIR / "alignments"
ANNOTATIONS_DIR = DEMO_DATA_DIR / "annotations"
CONFIG_DIR = DEMO_DATA_DIR / "config"
REFERENCE_DIR = DEMO_DATA_DIR / "reference"
SIGNALS_DIR = DEMO_DATA_DIR / "signals"
VARIANTS_DIR = DEMO_DATA_DIR / "variants"
REFERENCE_PATH = REFERENCE_DIR / "demo_reference.fa"
CNV_REFERENCE_PATH = REFERENCE_DIR / "demo_cnv_reference.fa"
LONG_REFERENCE_PATH = REFERENCE_DIR / "demo_long_reference.fa"
RNA_REFERENCE_PATH = REFERENCE_DIR / "demo_rna_fusion_reference.fa.gz"

CNV_CHROM = "chrCNV"
CNV_LENGTH = 80_000
CNV_READ_LENGTH = 150
CNV_BASELINE_DEPTH = 30.0
CNV_TUMOUR_PURITY = 0.75
SEQUENCING_ERROR_RATE = 0.0015
LONG_CHROM = "chrLong"
LONG_LENGTH = 12_000
RNA_CHROMS = ("chr2",)
# The reference only needs to extend beyond EML4 for fetching sequence, while
# BAM headers retain the GRCh37 chromosome-2 length so the ideogram is honest.
RNA_LENGTH = 42_530_000
GRCH37_CHR2_LENGTH = 243_199_373
EML4_DNA_BREAKPOINT = 42_523_383
ALK_DNA_BREAKPOINT = 29_448_041
EML4_EXONS = {
    12: (42_518_200, 42_518_320),
    13: (42_522_264, 42_522_399),
    14: (42_524_430, 42_524_560),
}
ALK_EXONS = {
    21: (29_445_100, 29_445_220),
    20: (29_446_205, 29_446_396),
    19: (29_448_324, 29_448_433),
}


@dataclass(frozen=True)
class CopyNumberSegment:
    start: int
    end: int
    label: str
    total_copy_number: int
    minor_copy_number: int

    def observed_copy_ratio(self, purity: float = CNV_TUMOUR_PURITY) -> float:
        tumour_copies = purity * self.total_copy_number
        normal_copies = (1.0 - purity) * 2
        return (tumour_copies + normal_copies) / 2

    def log2_ratio(self, purity: float = CNV_TUMOUR_PURITY) -> float:
        return log2(self.observed_copy_ratio(purity))

    def baf_centres(self, purity: float = CNV_TUMOUR_PURITY):
        total = purity * self.total_copy_number + (1.0 - purity) * 2
        minor = purity * self.minor_copy_number + (1.0 - purity)
        lower = minor / total
        if abs(lower - 0.5) < 1e-9:
            return (0.5,)
        return (lower, 1.0 - lower)


CNV_SEGMENTS = (
    CopyNumberSegment(0, 15_000, "CN2 diploid", 2, 1),
    CopyNumberSegment(15_000, 35_000, "CN1 loss + LOH", 1, 0),
    CopyNumberSegment(35_000, 50_000, "CN2 diploid", 2, 1),
    CopyNumberSegment(50_000, 70_000, "CN3 gain (1+2)", 3, 1),
    CopyNumberSegment(70_000, 80_000, "CN2 diploid", 2, 1),
)


@dataclass
class SimulatedVariant:
    position: int
    ref: str
    alt: str
    target_baf: float
    segment: CopyNumberSegment
    ref_depth: int = 0
    alt_depth: int = 0


def ensure_demo_directories() -> None:
    for directory in (
        ALIGNMENTS_DIR, ANNOTATIONS_DIR, CONFIG_DIR,
        REFERENCE_DIR, SIGNALS_DIR, VARIANTS_DIR,
    ):
        directory.mkdir(parents=True, exist_ok=True)


def alternate_base(reference_base: str) -> str:
    substitutions = {"A": "G", "C": "T", "G": "A", "T": "C"}
    return substitutions.get(reference_base, "A")


def add_sequencing_errors(sequence, qualities, rng: Random, error_rate: float) -> None:
    """Add sparse low-quality substitutions without changing read length."""
    for index, base in enumerate(sequence):
        if rng.random() >= error_rate:
            continue
        sequence[index] = rng.choice([candidate for candidate in "ACGT" if candidate != base])
        qualities[index] = 12


def create_variant_read(
    header, reference: str, read_index: int, sample_name: str,
    variant_fractions, rng: Random, read_group: str = None,
):
    read_length = 100
    start = max(45, min(100, round(rng.gauss(73, 11))))
    sequence = list(reference[start:start + read_length])
    for position, fraction in variant_fractions.items():
        query_position = position - start
        if query_position < 0 or query_position >= len(sequence):
            continue
        if rng.random() < fraction:
            sequence[query_position] = alternate_base(sequence[query_position])
    qualities = [40] * len(sequence)
    add_sequencing_errors(sequence, qualities, rng, SEQUENCING_ERROR_RATE)

    read = pysam.AlignedSegment(header)
    read.query_name = f"{sample_name}_read_{read_index + 1:03d}"
    read.query_sequence = "".join(sequence)
    read.flag = 16 if rng.random() < 0.5 else 0
    read.reference_id = 0
    read.reference_start = start
    read.mapping_quality = max(35, min(60, round(rng.gauss(56, 4))))
    read.cigar = [(0, read_length)]
    read.query_qualities = qualities
    if read_index % 3 != 2:
        read.set_tag("HP", 1 if read_index % 3 == 0 else 2)
        read.set_tag("PS", 1001 if read_index % 6 < 3 else 1002)
    if read_group:
        read.set_tag("RG", read_group)
    return read


def write_variant_bam(
    path: Path, sample_name: str, read_count: int, profile, seed: int,
    read_groups=None,
) -> None:
    with pysam.FastaFile(str(REFERENCE_PATH)) as fasta:
        reference = fasta.fetch("chrDemo").upper()
    rng = Random(seed)
    read_groups = list(read_groups or [])
    header_read_groups = (
        [{"ID": group, "SM": sample_name, "LB": f"{sample_name}_capture"}
         for group in read_groups]
        if read_groups else [{"ID": sample_name, "SM": sample_name}]
    )
    header_dict = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": "chrDemo", "LN": len(reference)}],
        "RG": header_read_groups,
    }
    header = pysam.AlignmentHeader.from_dict(header_dict)
    reads = []
    for read_index in range(read_count):
        reads.append(
            create_variant_read(
                header, reference, read_index, sample_name.lower(), profile, rng,
                read_group=(read_groups[read_index % len(read_groups)] if read_groups else None),
            )
        )
    reads.sort(key=lambda read: (read.reference_start, read.query_name))
    with pysam.AlignmentFile(str(path), "wb", header=header) as bam:
        for read in reads:
            bam.write(read)
    index_path = Path(f"{path}.bai")
    if index_path.exists():
        index_path.unlink()
    pysam.index(str(path))


def write_molecule_bam(path: Path, reference: str) -> None:
    """Simulate positional UMI families with PCR duplicates and sparse errors."""
    header = pysam.AlignmentHeader.from_dict({
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": "chrDemo", "LN": len(reference)}],
        "RG": [{"ID": "UMI", "SM": "UMI_consensus_demo"}],
    })
    rng = Random(20_026)
    variant_position = 118
    reads = []
    for family_index in range(26):
        umi = f"UMI{family_index + 1:04d}"
        cell = "CELL-A" if family_index < 18 else "CELL-B"
        family_size = 2 + (family_index * 5) % 6
        family_start = 48 + (family_index * 7) % 34
        carries_variant = family_index in {1, 3, 4, 8, 11, 15, 18, 21, 24}
        for copy_index in range(family_size):
            start = family_start + (1 if copy_index == family_size - 1 and family_index % 4 == 0 else 0)
            sequence = list(reference[start:start + 100])
            if carries_variant:
                query_position = variant_position - start
                sequence[query_position] = alternate_base(reference[variant_position])
            qualities = [38] * len(sequence)
            add_sequencing_errors(sequence, qualities, rng, 0.0025)

            read = pysam.AlignedSegment(header)
            read.query_name = f"{umi}_copy_{copy_index + 1:02d}"
            read.query_sequence = "".join(sequence)
            reverse = family_index % 5 == 0 and copy_index == family_size - 1
            read.flag = (16 if reverse else 0) | (1024 if copy_index else 0)
            read.reference_id = 0
            read.reference_start = start
            read.mapping_quality = max(42, min(60, round(rng.gauss(57, 3))))
            read.cigar = [(0, len(sequence))]
            read.query_qualities = qualities
            read.set_tag("RG", "UMI")
            read.set_tag("RX", umi)
            read.set_tag("CB", cell)
            reads.append(read)

    # A few true singleton molecules illustrate optional family-size filtering.
    for singleton_index in range(3):
        start = 54 + singleton_index * 9
        sequence = reference[start:start + 100]
        read = pysam.AlignedSegment(header)
        read.query_name = f"singleton_{singleton_index + 1}"
        read.query_sequence = sequence
        read.flag = 0
        read.reference_id = 0
        read.reference_start = start
        read.mapping_quality = 55
        read.cigar = [(0, len(sequence))]
        read.query_qualities = pysam.qualitystring_to_array("I" * len(sequence))
        read.set_tag("RG", "UMI")
        read.set_tag("RX", f"SINGLE{singleton_index + 1}")
        read.set_tag("CB", "CELL-B")
        reads.append(read)

    reads.sort(key=lambda read: (read.reference_start, read.query_name))
    with pysam.AlignmentFile(str(path), "wb", header=header) as bam:
        for read in reads:
            bam.write(read)
    index_path = Path(f"{path}.bai")
    if index_path.exists():
        index_path.unlink()
    pysam.index(str(path))


def write_insertion_bam(path: Path) -> None:
    """Write close-zoom reads with short CIGAR insertions at one breakpoint."""
    with pysam.FastaFile(str(REFERENCE_PATH)) as fasta:
        reference = fasta.fetch("chrDemo").upper()
    header_dict = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": "chrDemo", "LN": len(reference)}],
        "RG": [{"ID": "INS", "SM": "Insertion_demo"}],
    }
    header = pysam.AlignmentHeader.from_dict(header_dict)
    reads = []
    breakpoint = 120
    motifs = ("TGCA", "ACGT", "GGTT", "CAGA", "TTAC", "AGGC")
    for read_index, inserted_sequence in enumerate(motifs):
        start = 91 + read_index
        left_match = breakpoint - start
        right_match = 55
        sequence = (
            reference[start:breakpoint]
            + inserted_sequence
            + reference[breakpoint:breakpoint + right_match]
        )
        read = pysam.AlignedSegment(header)
        read.query_name = f"insertion_read_{read_index + 1:02d}"
        read.query_sequence = sequence
        read.flag = 16 if read_index % 2 else 0
        read.reference_id = 0
        read.reference_start = start
        read.mapping_quality = 60
        read.cigar = [(0, left_match), (1, len(inserted_sequence)), (0, right_match)]
        read.query_qualities = pysam.qualitystring_to_array("I" * len(sequence))
        read.set_tag("RG", "INS")
        reads.append(read)
    reads.sort(key=lambda read: (read.reference_start, read.query_name))
    with pysam.AlignmentFile(str(path), "wb", header=header) as bam:
        for read in reads:
            bam.write(read)
    index_path = Path(f"{path}.bai")
    if index_path.exists():
        index_path.unlink()
    pysam.index(str(path))


def copy_number_segment_at(position: int) -> CopyNumberSegment:
    for segment in CNV_SEGMENTS:
        if segment.start <= position < segment.end:
            return segment
    raise ValueError(f"Position {position} is outside the simulated CNV locus.")


def write_cnv_reference(path: Path) -> str:
    rng = Random(18_018)
    sequence = "".join(rng.choices("ACGT", weights=(29, 21, 21, 29), k=CNV_LENGTH))
    with path.open("w", encoding="utf-8") as handle:
        handle.write(f">{CNV_CHROM}\n")
        for offset in range(0, len(sequence), 80):
            handle.write(sequence[offset:offset + 80] + "\n")
    index_path = Path(f"{path}.fai")
    if index_path.exists():
        index_path.unlink()
    pysam.faidx(str(path))
    return sequence


def write_small_reference(path: Path) -> str:
    rng = Random(8_008)
    sequence = "".join(rng.choices("ACGT", weights=(29, 21, 21, 29), k=800))
    with path.open("w", encoding="utf-8") as handle:
        handle.write(">chrDemo\n")
        for offset in range(0, len(sequence), 80):
            handle.write(sequence[offset:offset + 80] + "\n")
    index_path = Path(f"{path}.fai")
    if index_path.exists():
        index_path.unlink()
    pysam.faidx(str(path))
    return sequence


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def write_long_reference(path: Path) -> str:
    """Write a GC-balanced synthetic locus for ONT-style long reads."""
    rng = Random(88_021)
    sequence = "".join(rng.choices("ACGT", weights=(26, 24, 24, 26), k=LONG_LENGTH))
    with path.open("w", encoding="utf-8") as handle:
        handle.write(f">{LONG_CHROM}\n")
        for offset in range(0, len(sequence), 80):
            handle.write(sequence[offset:offset + 80] + "\n")
    index_path = Path(f"{path}.fai")
    if index_path.exists():
        index_path.unlink()
    pysam.faidx(str(path))
    return sequence


def write_rna_fusion_reference(path: Path) -> dict[str, str]:
    """Write a compact GRCh37-like chr2 reference for the EML4::ALK demo."""
    pattern = b"ACGTTGCAAGTCGATCGGATCCATGCTA"
    sequence = bytearray((pattern * (RNA_LENGTH // len(pattern) + 1))[:RNA_LENGTH])

    # Make the visible windows locally non-repetitive without constructing a
    # chromosome-sized list of Python strings.
    rng = Random(24_003)
    for start, end in (
        (ALK_EXONS[21][0] - 1_000, ALK_EXONS[19][1] + 1_000),
        (EML4_EXONS[12][0] - 1_000, EML4_EXONS[14][1] + 1_000),
    ):
        sequence[start:end] = bytes(
            ord(base) for base in rng.choices("ACGT", weights=(29, 21, 21, 29), k=end - start)
        )

    # Canonical splice motifs. EML4 is on +; ALK is on the reverse strand, so
    # a transcript-oriented GT-AG intron appears as genomic CT...AC.
    sequence[EML4_EXONS[12][1]:EML4_EXONS[12][1] + 2] = b"GT"
    sequence[EML4_EXONS[13][0] - 2:EML4_EXONS[13][0]] = b"AG"
    sequence[EML4_EXONS[13][1]:EML4_EXONS[13][1] + 2] = b"GT"
    sequence[EML4_EXONS[14][0] - 2:EML4_EXONS[14][0]] = b"AG"
    for lower_exon, upper_exon in ((21, 20), (20, 19)):
        sequence[ALK_EXONS[lower_exon][1]:ALK_EXONS[lower_exon][1] + 2] = b"CT"
        sequence[ALK_EXONS[upper_exon][0] - 2:ALK_EXONS[upper_exon][0]] = b"AC"

    reference = sequence.decode("ascii")
    with pysam.BGZFile(str(path), "w") as handle:
        handle.write(b">chr2\n")
        for offset in range(0, len(reference), 80):
            handle.write((reference[offset:offset + 80] + "\n").encode("ascii"))
    for suffix in (".fai", ".gzi"):
        index_path = Path(f"{path}{suffix}")
        if index_path.exists():
            index_path.unlink()
    pysam.faidx(str(path))
    return {"chr2": reference}


def write_rna_fusion_gtf(path: Path) -> None:
    genes = (
        ("+", "EML4", "EML4-201", EML4_EXONS),
        ("-", "ALK", "ALK-201", ALK_EXONS),
    )
    with path.open("w", encoding="utf-8") as handle:
        for strand, gene, transcript, numbered_exons in genes:
            exons = sorted(numbered_exons.values())
            attributes = f'gene_id "{gene}"; gene_name "{gene}";'
            handle.write(
                f"chr2\tdemo\tgene\t{exons[0][0] + 1}\t{exons[-1][1]}\t.\t"
                f"{strand}\t.\t{attributes}\n"
            )
            transcript_attributes = (
                attributes + f' transcript_id "{transcript}"; '
                f'transcript_name "{transcript}";'
            )
            handle.write(
                f"chr2\tdemo\ttranscript\t{exons[0][0] + 1}\t{exons[-1][1]}\t.\t"
                f"{strand}\t.\t{transcript_attributes}\n"
            )
            for exon_number, (start, end) in numbered_exons.items():
                exon_attributes = transcript_attributes + f' exon_number "{exon_number}";'
                handle.write(
                    f"chr2\tdemo\texon\t{start + 1}\t{end}\t.\t{strand}\t.\t"
                    f"{exon_attributes}\n"
                )


def write_rna_fusion_bam(path: Path, references: dict[str, str]) -> None:
    """Simulate EML4 exon 13::ALK exon 20 (variant 1) RNA evidence."""
    reference = references["chr2"]
    header = pysam.AlignmentHeader.from_dict({
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": "chr2", "LN": GRCH37_CHR2_LENGTH}],
        "RG": [{"ID": "RNA", "SM": "EML4_ALK_V1_RNA", "PL": "ILLUMINA"}],
    })
    rng = Random(24_004)
    reads = []

    def sequence_with_errors(start, end):
        sequence = list(reference[start:end])
        qualities = [38] * len(sequence)
        add_sequencing_errors(sequence, qualities, rng, SEQUENCING_ERROR_RATE)
        return "".join(sequence), qualities

    def add_spliced(name, left_start, donor, acceptor, count, reverse=False):
        left_length = donor - left_start
        right_length = 60
        intron = acceptor - donor
        for index in range(count):
            jitter = index % 5
            start = left_start - jitter
            adjusted_left = left_length + jitter
            sequence = list(reference[start:donor] + reference[acceptor:acceptor + right_length])
            qualities = [38] * len(sequence)
            add_sequencing_errors(sequence, qualities, rng, SEQUENCING_ERROR_RATE)
            read = pysam.AlignedSegment(header)
            read.query_name = f"{name}_{index + 1:03d}"
            read.query_sequence = "".join(sequence)
            read.flag = 16 if reverse else 0
            read.reference_id = 0
            read.reference_start = start
            read.mapping_quality = max(48, min(60, round(rng.gauss(58, 2))))
            read.cigar = [(0, adjusted_left), (3, intron), (0, right_length)]
            read.query_qualities = qualities
            read.set_tag("RG", "RNA")
            read.set_tag("NH", 1)
            reads.append(read)

    add_spliced("EML4_e12_e13", EML4_EXONS[12][1] - 60,
                EML4_EXONS[12][1], EML4_EXONS[13][0], 34)
    add_spliced("EML4_e13_e14", EML4_EXONS[13][1] - 55,
                EML4_EXONS[13][1], EML4_EXONS[14][0], 8)
    add_spliced("ALK_e21_e20", ALK_EXONS[21][1] - 60,
                ALK_EXONS[21][1], ALK_EXONS[20][0], 24, reverse=True)
    add_spliced("ALK_e20_e19", ALK_EXONS[20][1] - 55,
                ALK_EXONS[20][1], ALK_EXONS[19][0], 28, reverse=True)

    # Eighteen split fragments cross the mature V1 junction. The EML4 segment
    # ends at exon 13 and the reverse-strand ALK segment ends at exon 20.
    eml4_rna_breakpoint = EML4_EXONS[13][1]
    alk_rna_breakpoint = ALK_EXONS[20][1]
    for index in range(18):
        jitter = index % 5 - 2
        first_start = eml4_rna_breakpoint - 55 + jitter
        second_start = alk_rna_breakpoint - 55 + jitter
        first_sequence = reference[first_start:first_start + 55]
        # BAM stores aligned bases in reference orientation even when FLAG 0x10
        # records the segment's reverse strand.
        second_sequence = reference[second_start:second_start + 55]
        query_sequence = list(first_sequence + second_sequence)
        qualities = [38] * len(query_sequence)
        add_sequencing_errors(query_sequence, qualities, rng, SEQUENCING_ERROR_RATE)
        query_sequence = "".join(query_sequence)
        supplementary_sequence = list(second_sequence + first_sequence)
        supplementary_qualities = [38] * len(supplementary_sequence)
        add_sequencing_errors(
            supplementary_sequence, supplementary_qualities, rng,
            SEQUENCING_ERROR_RATE,
        )
        supplementary_sequence = "".join(supplementary_sequence)
        query_name = f"EML4e13_ALKe20_split_{index + 1:02d}"

        primary = pysam.AlignedSegment(header)
        primary.query_name = query_name
        primary.query_sequence = query_sequence
        primary.flag = 0
        primary.reference_id = 0
        primary.reference_start = first_start
        primary.mapping_quality = 60
        primary.cigar = [(0, 55), (4, 55)]
        primary.query_qualities = qualities
        primary.set_tag("SA", f"chr2,{second_start + 1},-,55M55S,58,1;")
        primary.set_tag("RG", "RNA")
        primary.set_tag("NH", 1)
        reads.append(primary)

        supplementary = pysam.AlignedSegment(header)
        supplementary.query_name = query_name
        supplementary.query_sequence = supplementary_sequence
        supplementary.flag = 2048 | 16
        supplementary.reference_id = 0
        supplementary.reference_start = second_start
        supplementary.mapping_quality = 58
        supplementary.cigar = [(0, 55), (4, 55)]
        supplementary.query_qualities = supplementary_qualities
        supplementary.set_tag("SA", f"chr2,{first_start + 1},+,55M55S,60,1;")
        supplementary.set_tag("RG", "RNA")
        supplementary.set_tag("NH", 1)
        reads.append(supplementary)

    # Ten spanning fragments add paired-end support without either read itself
    # crossing the chimeric junction.
    for index in range(10):
        first_start = eml4_rna_breakpoint - 100 - index % 4
        second_start = alk_rna_breakpoint - 100 - index % 4
        query_name = f"EML4e13_ALKe20_pair_{index + 1:02d}"
        for read_side, start, mate_start, reverse, mate_reverse in (
            (1, first_start, second_start, False, False),
            (2, second_start, first_start, False, False),
        ):
            sequence, qualities = sequence_with_errors(start, start + 100)
            read = pysam.AlignedSegment(header)
            read.query_name = query_name
            read.query_sequence = sequence
            read.flag = (
                1 | (64 if read_side == 1 else 128)
                | (16 if reverse else 0)
                | (32 if mate_reverse else 0)
            )
            read.reference_id = 0
            read.reference_start = start
            read.mapping_quality = 55
            read.cigar = [(0, 100)]
            read.next_reference_id = 0
            read.next_reference_start = mate_start
            read.template_length = 0
            read.query_qualities = qualities
            read.set_tag("RG", "RNA")
            read.set_tag("NH", 1)
            reads.append(read)

    reads.sort(key=lambda read: (read.reference_id, read.reference_start, read.query_name))
    with pysam.AlignmentFile(str(path), "wb", header=header) as bam:
        for read in reads:
            bam.write(read)
    index_path = Path(f"{path}.bai")
    if index_path.exists():
        index_path.unlink()
    pysam.index(str(path))


def write_eml4_alk_dna_bam(path: Path, references: dict[str, str]) -> None:
    """Simulate heterozygous DNA evidence for the EML4::ALK inversion."""
    reference = references["chr2"]
    header = pysam.AlignmentHeader.from_dict({
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": "chr2", "LN": GRCH37_CHR2_LENGTH}],
        "RG": [{"ID": "DNA", "SM": "EML4_ALK_inversion_DNA", "PL": "ILLUMINA"}],
    })
    rng = Random(24_005)
    reads = []

    def make_read(name, start, length=150, reverse=False, mapq=60):
        sequence = list(reference[start:start + length])
        qualities = [38] * len(sequence)
        add_sequencing_errors(sequence, qualities, rng, SEQUENCING_ERROR_RATE)
        read = pysam.AlignedSegment(header)
        read.query_name = name
        read.query_sequence = "".join(sequence)
        read.flag = 16 if reverse else 0
        read.reference_id = 0
        read.reference_start = start
        read.mapping_quality = mapq
        read.cigar = [(0, length)]
        read.query_qualities = qualities
        read.set_tag("RG", "DNA")
        return read

    # Approximately 15x local reference support from ordinary 350-bp FR
    # fragments, with start and MAPQ variation and a sparse error model.
    for locus_name, window_start, window_end in (
        ("ALK", ALK_DNA_BREAKPOINT - 1_650, ALK_DNA_BREAKPOINT + 1_650),
        ("EML4", EML4_DNA_BREAKPOINT - 1_450, EML4_DNA_BREAKPOINT + 1_450),
    ):
        for index in range(132):
            fragment_length = max(300, min(430, round(rng.gauss(355, 30))))
            first_start = rng.randint(window_start, window_end - fragment_length)
            second_start = first_start + fragment_length - 150
            query_name = f"DNA_{locus_name}_concordant_{index + 1:03d}"
            first = make_read(query_name, first_start, mapq=rng.randint(48, 60))
            second = make_read(query_name, second_start, reverse=True, mapq=rng.randint(48, 60))
            first.flag = 1 | 2 | 64 | 32
            second.flag = 1 | 2 | 128 | 16
            first.next_reference_id = second.next_reference_id = 0
            first.next_reference_start = second_start
            second.next_reference_start = first_start
            first.template_length = fragment_length
            second.template_length = -fragment_length
            reads.extend((first, second))

    # Split reads directly resolve the two intronic DNA breakpoints. Their
    # opposite SA strands encode the inversion rather than a deletion.
    for index in range(10):
        jitter = index % 5 - 2
        eml4_start = EML4_DNA_BREAKPOINT - 60 + jitter
        alk_start = ALK_DNA_BREAKPOINT + jitter
        eml4_sequence = reference[eml4_start:eml4_start + 60]
        alk_sequence = reference[alk_start:alk_start + 60]
        query = list(eml4_sequence + alk_sequence)
        qualities = [38] * len(query)
        add_sequencing_errors(query, qualities, rng, SEQUENCING_ERROR_RATE)
        supplementary_query = list(alk_sequence + eml4_sequence)
        supplementary_qualities = [38] * len(supplementary_query)
        add_sequencing_errors(
            supplementary_query, supplementary_qualities, rng,
            SEQUENCING_ERROR_RATE,
        )
        query_name = f"DNA_inv2_split_{index + 1:02d}"

        primary = pysam.AlignedSegment(header)
        primary.query_name = query_name
        primary.query_sequence = "".join(query)
        primary.flag = 0
        primary.reference_id = 0
        primary.reference_start = eml4_start
        primary.mapping_quality = 60
        primary.cigar = [(0, 60), (4, 60)]
        primary.query_qualities = qualities
        primary.set_tag("SA", f"chr2,{alk_start + 1},-,60S60M,58,1;")
        primary.set_tag("RG", "DNA")
        reads.append(primary)

        supplementary = pysam.AlignedSegment(header)
        supplementary.query_name = query_name
        supplementary.query_sequence = "".join(supplementary_query)
        supplementary.flag = 2048 | 16
        supplementary.reference_id = 0
        supplementary.reference_start = alk_start
        supplementary.mapping_quality = 58
        supplementary.cigar = [(4, 60), (0, 60)]
        supplementary.query_qualities = supplementary_qualities
        supplementary.set_tag("SA", f"chr2,{eml4_start + 1},+,60M60S,60,1;")
        supplementary.set_tag("RG", "DNA")
        reads.append(supplementary)

    # Same-strand pairs are the characteristic short-read geometry of an
    # inversion. Model both reciprocal junction orientations (FF and RR).
    for orientation, count in (("FF", 10), ("RR", 7)):
        for index in range(count):
            jitter = index % 5 - 2
            if orientation == "FF":
                eml4_start = EML4_DNA_BREAKPOINT - 150 + jitter
                alk_start = ALK_DNA_BREAKPOINT - 150 + jitter
                reverse = False
            else:
                eml4_start = EML4_DNA_BREAKPOINT + jitter
                alk_start = ALK_DNA_BREAKPOINT + jitter
                reverse = True
            query_name = f"DNA_inv2_{orientation}_{index + 1:02d}"
            eml4 = make_read(query_name, eml4_start, reverse=reverse, mapq=56)
            alk = make_read(query_name, alk_start, reverse=reverse, mapq=56)
            eml4.flag = 1 | 64 | (16 if reverse else 0) | (32 if reverse else 0)
            alk.flag = 1 | 128 | (16 if reverse else 0) | (32 if reverse else 0)
            eml4.next_reference_id = alk.next_reference_id = 0
            eml4.next_reference_start = alk_start
            alk.next_reference_start = eml4_start
            eml4.template_length = alk.template_length = 0
            reads.extend((eml4, alk))

    reads.sort(key=lambda read: (read.reference_id, read.reference_start, read.query_name, read.flag))
    with pysam.AlignmentFile(str(path), "wb", header=header) as bam:
        for read in reads:
            bam.write(read)
    index_path = Path(f"{path}.bai")
    if index_path.exists():
        index_path.unlink()
    pysam.index(str(path))


def encode_modification_deltas(sequence: str, canonical_base: str, selected) -> list[int]:
    """Encode selected query indices as SAM MM canonical-base deltas."""
    canonical_positions = [
        index for index, base in enumerate(sequence) if base == canonical_base
    ]
    ordinal_by_position = {
        position: ordinal for ordinal, position in enumerate(canonical_positions)
    }
    deltas = []
    previous_ordinal = -1
    for position in sorted(selected):
        ordinal = ordinal_by_position[position]
        deltas.append(ordinal - previous_ordinal - 1)
        previous_ordinal = ordinal
    return deltas


def write_long_read_bam(path: Path, reference: str) -> None:
    """Simulate ONT-like reads with indels, sparse errors, and MM/ML calls."""
    header = pysam.AlignmentHeader.from_dict({
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": LONG_CHROM, "LN": len(reference)}],
        "RG": [{"ID": "ONT", "SM": "Long_read_modifications", "PL": "ONT"}],
    })
    rng = Random(88_022)
    reads = []
    for read_index in range(26):
        start = 2200 + rng.randint(0, 1900)
        aligned_span = rng.randint(3600, 5200)
        end = min(start + aligned_span, len(reference) - 50)
        aligned_span = end - start
        event_position = rng.randint(900, max(901, aligned_span - 900))
        event_length = rng.randint(3, 18)
        if read_index % 3 == 0:
            left = reference[start:start + event_position]
            inserted = "".join(rng.choices("ACGT", k=event_length))
            right = reference[start + event_position:end]
            sequence = list(left + inserted + right)
            cigar = [(0, event_position), (1, event_length), (0, aligned_span - event_position)]
        elif read_index % 3 == 1:
            deletion_length = min(event_length, aligned_span - event_position - 1)
            left = reference[start:start + event_position]
            right = reference[start + event_position + deletion_length:end]
            sequence = list(left + right)
            cigar = [
                (0, event_position), (2, deletion_length),
                (0, aligned_span - event_position - deletion_length),
            ]
        else:
            sequence = list(reference[start:end])
            cigar = [(0, aligned_span)]

        qualities = [max(18, min(40, round(rng.gauss(31, 4)))) for _ in sequence]
        add_sequencing_errors(sequence, qualities, rng, 0.0035)
        stored_sequence = "".join(sequence)
        reverse = read_index % 4 == 1
        original_sequence = reverse_complement(stored_sequence) if reverse else stored_sequence

        c_candidates = [
            index for index in range(len(original_sequence) - 1)
            if original_sequence[index:index + 2] == "CG"
        ]
        a_candidates = [
            index for index, base in enumerate(original_sequence) if base == "A"
        ]
        c_selected = [
            position for position in c_candidates
            if rng.random() < (0.42 if read_index < 14 else 0.18)
        ]
        a_selected = [position for position in a_candidates if rng.random() < 0.012]
        mm_groups = []
        ml_values = []
        if c_selected:
            deltas = encode_modification_deltas(original_sequence, "C", c_selected)
            mm_groups.append("C+m.," + ",".join(str(value) for value in deltas) + ";")
            ml_values.extend(rng.randint(185, 255) for _ in c_selected)
        if a_selected:
            deltas = encode_modification_deltas(original_sequence, "A", a_selected)
            mm_groups.append("A+a.," + ",".join(str(value) for value in deltas) + ";")
            ml_values.extend(rng.randint(165, 245) for _ in a_selected)

        read = pysam.AlignedSegment(header)
        read.query_name = f"ont_read_{read_index + 1:03d}"
        read.query_sequence = stored_sequence
        read.flag = (16 if reverse else 0) | (2048 if read_index in (6, 19) else 0)
        read.reference_id = 0
        read.reference_start = start
        read.mapping_quality = max(20, min(60, round(rng.gauss(52, 7))))
        read.cigar = cigar
        read.query_qualities = qualities
        read.set_tag("RG", "ONT")
        if mm_groups:
            read.set_tag("MM", "".join(mm_groups))
            read.set_tag("ML", array("B", ml_values))
        if read.is_supplementary:
            read.set_tag("SA", f"{LONG_CHROM},{start + 2401},+,1200M,40,8;")
        reads.append(read)

    reads.sort(key=lambda read: (read.reference_start, read.query_name))
    with pysam.AlignmentFile(str(path), "wb", header=header) as bam:
        for read in reads:
            bam.write(read)
    index_path = Path(f"{path}.bai")
    if index_path.exists():
        index_path.unlink()
    pysam.index(str(path))


def write_small_variant_vcf(
    path: Path, reference: str, profile, sample_prefix: str,
) -> None:
    with path.open("w", encoding="utf-8") as handle:
        handle.write("##fileformat=VCFv4.2\n")
        handle.write(f"##contig=<ID=chrDemo,length={len(reference)}>\n")
        handle.write(
            '##INFO=<ID=EXPECTED_AF,Number=1,Type=Float,'
            'Description="Expected synthetic allele fraction">\n'
        )
        handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for index, (position, fraction) in enumerate(sorted(profile.items()), start=1):
            ref = reference[position]
            handle.write(
                f"chrDemo\t{position + 1}\t{sample_prefix}-SNV-{index:02d}\t"
                f"{ref}\t{alternate_base(ref)}\t100\tPASS\t"
                f"EXPECTED_AF={fraction:.3f}\n"
            )


def build_cnv_variant_sites(reference: str) -> list[SimulatedVariant]:
    rng = Random(18_019)
    sites = []
    cursor = 700
    while cursor < CNV_LENGTH - 500:
        position = min(cursor + rng.randint(-160, 160), CNV_LENGTH - 501)
        segment = copy_number_segment_at(position)
        centres = segment.baf_centres()
        target_baf = centres[len(sites) % len(centres)]
        ref = reference[position]
        sites.append(SimulatedVariant(
            position=position,
            ref=ref,
            alt=alternate_base(ref),
            target_baf=target_baf,
            segment=segment,
        ))
        cursor += rng.randint(780, 1_120)
    return sites


def write_cnv_tumour_bam(
    path: Path, reference: str, sites: list[SimulatedVariant],
) -> int:
    """Simulate one impure tumour whose depth and BAF follow the CN state."""
    header_dict = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": CNV_CHROM, "LN": len(reference)}],
        "RG": [{"ID": "TumourCNV", "SM": "Tumour_CNV_75pct"}],
    }
    header = pysam.AlignmentHeader.from_dict(header_dict)
    rng = Random(18_020)
    reads = []
    read_index = 0
    for segment in CNV_SEGMENTS:
        target_depth = CNV_BASELINE_DEPTH * segment.observed_copy_ratio()
        read_count = round(
            (segment.end - segment.start) * target_depth / CNV_READ_LENGTH
        )
        segment_sites = [
            site for site in sites
            if segment.start <= site.position < segment.end
        ]
        for _ in range(read_count):
            read_index += 1
            midpoint = rng.randrange(segment.start, segment.end)
            start = max(0, min(
                len(reference) - CNV_READ_LENGTH,
                midpoint - CNV_READ_LENGTH // 2,
            ))
            sequence = list(reference[start:start + CNV_READ_LENGTH])
            for site in segment_sites:
                query_position = site.position - start
                if 0 <= query_position < len(sequence) and rng.random() < site.target_baf:
                    sequence[query_position] = site.alt
            qualities = [36] * len(sequence)
            add_sequencing_errors(
                sequence, qualities, rng, SEQUENCING_ERROR_RATE
            )
            for site in segment_sites:
                query_position = site.position - start
                if not 0 <= query_position < len(sequence):
                    continue
                if sequence[query_position] == site.ref:
                    site.ref_depth += 1
                elif sequence[query_position] == site.alt:
                    site.alt_depth += 1

            read = pysam.AlignedSegment(header)
            read.query_name = f"tumour_cnv_read_{read_index:05d}"
            read.query_sequence = "".join(sequence)
            read.flag = 16 if rng.random() < 0.5 else 0
            read.reference_id = 0
            read.reference_start = start
            read.mapping_quality = max(35, min(60, round(rng.gauss(55, 5))))
            read.cigar = [(0, CNV_READ_LENGTH)]
            read.query_qualities = qualities
            read.set_tag("RG", "TumourCNV")
            reads.append(read)

    reads.sort(key=lambda read: (read.reference_start, read.query_name))
    with pysam.AlignmentFile(str(path), "wb", header=header) as bam:
        for read in reads:
            bam.write(read)
    index_path = Path(f"{path}.bai")
    if index_path.exists():
        index_path.unlink()
    pysam.index(str(path))
    return len(reads)


def write_cnv_baf_vcf(path: Path, sites: list[SimulatedVariant]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        handle.write("##fileformat=VCFv4.2\n")
        handle.write(f"##contig=<ID={CNV_CHROM},length={CNV_LENGTH}>\n")
        handle.write('##INFO=<ID=STATE,Number=1,Type=String,Description="Copy-number state">\n')
        handle.write('##INFO=<ID=TCN,Number=1,Type=Integer,Description="Tumour total copy number">\n')
        handle.write('##INFO=<ID=MCN,Number=1,Type=Integer,Description="Tumour minor copy number">\n')
        handle.write('##INFO=<ID=TARGET_BAF,Number=1,Type=Float,Description="Purity-adjusted expected BAF">\n')
        handle.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        handle.write('##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allele depths">\n')
        handle.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">\n')
        handle.write(
            f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTumour\n"
        )
        for index, site in enumerate(sites, start=1):
            depth = site.ref_depth + site.alt_depth
            state = site.segment.label.replace(" ", "_").replace("+", "plus")
            handle.write(
                f"{CNV_CHROM}\t{site.position + 1}\tcnv-het-{index:03d}\t"
                f"{site.ref}\t{site.alt}\t100\tPASS\tSTATE={state};"
                f"TCN={site.segment.total_copy_number};"
                f"MCN={site.segment.minor_copy_number};"
                f"TARGET_BAF={site.target_baf:.4f}\tGT:AD:DP\t"
                f"0/1:{site.ref_depth},{site.alt_depth}:{depth}\n"
            )


def write_copy_number_tracks(seg_path: Path, state_path: Path) -> None:
    legacy_segments = (
        (101_867_481, 101_867_500, 18, -0.70),
        (101_867_501, 101_867_520, 21, -0.42),
        (101_867_521, 101_867_540, 17, -0.08),
        (101_867_541, 101_867_560, 20, 0.28),
        (101_867_561, 101_867_580, 24, 0.62),
        (101_867_581, 101_867_600, 19, 0.88),
        (101_867_601, 101_867_620, 16, 0.18),
    )
    with seg_path.open("w", encoding="utf-8") as handle:
        handle.write("Sample\tChromosome\tStart\tEnd\tNum_Probes\tSegment_Mean\n")
        for start, end, probes, value in legacy_segments:
            handle.write(
                f"Tumour\tchr9\t{start}\t{end}\t{probes}\t{value:.2f}\n"
            )
        for segment in CNV_SEGMENTS:
            probe_count = max(20, (segment.end - segment.start) // 250)
            handle.write(
                f"Tumour\t{CNV_CHROM}\t{segment.start + 1}\t{segment.end}\t"
                f"{probe_count}\t{segment.log2_ratio():.3f}\n"
            )
    with state_path.open("w", encoding="utf-8") as handle:
        for segment in CNV_SEGMENTS:
            baf = "/".join(f"{centre:.2f}" for centre in segment.baf_centres())
            label = f"{segment.label} · BAF {baf}"
            handle.write(
                f"{CNV_CHROM}\t{segment.start}\t{segment.end}\t{label}\n"
            )


def create_rna_read(header, name: str, start: int, cigar, reverse: bool):
    query_length = 0
    for operation, length in cigar:
        if operation in (0, 1, 4, 7, 8):
            query_length += length
    read = pysam.AlignedSegment(header)
    read.query_name = name
    read.query_sequence = "A" * query_length
    read.flag = 16 if reverse else 0
    read.reference_id = 0
    read.reference_start = start
    read.mapping_quality = 60
    read.cigar = cigar
    read.query_qualities = pysam.qualitystring_to_array("I" * query_length)
    read.set_tag("XS", "-" if reverse else "+")
    return read


def write_met_ex14_bam(path: Path) -> None:
    """Write a synthetic METex14-positive RNA-seq cohort on GRCh38.

    Coordinates follow the MANE Select transcript NM_000245.4.  The dominant
    exon 13-to-15 junction models exon 14 skipping caused by c.3028+1G>T;
    lower-support exon-inclusion junctions remain as realistic background.
    """
    header_dict = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": "chr7", "LN": 159_345_973}],
        "RG": [{"ID": "METex14", "SM": "METex14_LUAD"}],
    }
    header = pysam.AlignmentHeader.from_dict(header_dict)
    reads = []

    # Zero-based, half-open equivalents of the GRCh38 exon coordinates:
    # exon 13 116771498-116771654, exon 14 116771849-116771989,
    # exon 15 116774881-116775111 (Ensembl canonical/MANE transcript).
    exon_13 = (116_771_497, 116_771_654)
    exon_14 = (116_771_848, 116_771_989)
    exon_15 = (116_774_880, 116_775_111)

    def add_exonic_reads(exon, count: int, label: str) -> None:
        read_length = 90
        available_starts = exon[1] - exon[0] - read_length + 1
        for read_index in range(count):
            start = exon[0] + (read_index * 17) % available_starts
            reads.append(create_rna_read(
                header, f"{label}_{read_index + 1:03d}", start,
                [(0, read_length)], False,
            ))

    def add_junction_reads(
        left_exon_end: int, right_exon_start: int, count: int, label: str,
    ) -> None:
        intron_length = right_exon_start - left_exon_end
        for read_index in range(count):
            left_match = 42 + (read_index * 7) % 17
            right_match = 100 - left_match
            reads.append(create_rna_read(
                header, f"{label}_{read_index + 1:03d}",
                left_exon_end - left_match,
                [(0, left_match), (3, intron_length), (0, right_match)], False,
            ))

    add_exonic_reads(exon_13, 56, "met_exon13")
    add_exonic_reads(exon_14, 24, "met_exon14")
    add_exonic_reads(exon_15, 72, "met_exon15")
    add_junction_reads(exon_13[1], exon_15[0], 96, "met_ex14_skipping")
    add_junction_reads(exon_13[1], exon_14[0], 28, "met_ex13_14_inclusion")
    add_junction_reads(exon_14[1], exon_15[0], 24, "met_ex14_15_inclusion")

    reads.sort(key=lambda read: (read.reference_start, read.query_name))
    with pysam.AlignmentFile(str(path), "wb", header=header) as bam:
        for read in reads:
            bam.write(read)
    index_path = Path(f"{path}.bai")
    if index_path.exists():
        index_path.unlink()
    pysam.index(str(path))


def write_met_ex14_tracks() -> None:
    """Write MET gene context and the canonical exon-14 donor SNV."""
    attributes = (
        'gene_id "MET"; transcript_id "NM_000245.4"; gene_name "MET"; '
        'tag "MANE_Select";'
    )
    with (ANNOTATIONS_DIR / "demo_met_ex14.gtf").open("w", encoding="utf-8") as handle:
        handle.write(
            f"chr7\tdemo\ttranscript\t116771498\t116775111\t.\t+\t.\t{attributes}\n"
        )
        for exon_number, (start, end) in enumerate([
            (116_771_498, 116_771_654),
            (116_771_849, 116_771_989),
            (116_774_881, 116_775_111),
        ], start=13):
            exon_attributes = f'{attributes} exon_number "{exon_number}";'
            handle.write(
                f"chr7\tdemo\texon\t{start}\t{end}\t.\t+\t.\t{exon_attributes}\n"
            )
            handle.write(
                f"chr7\tdemo\tCDS\t{start}\t{end}\t.\t+\t0\t{exon_attributes}\n"
            )

    with (VARIANTS_DIR / "demo_met_ex14.vcf").open("w", encoding="utf-8") as handle:
        handle.write("##fileformat=VCFv4.2\n")
        handle.write("##reference=GRCh38\n")
        handle.write("##contig=<ID=chr7,length=159345973>\n")
        handle.write(
            '##INFO=<ID=GENE,Number=1,Type=String,Description="Gene symbol">\n'
        )
        handle.write(
            '##INFO=<ID=RNA_EFFECT,Number=1,Type=String,Description="Observed RNA effect">\n'
        )
        handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        handle.write(
            "chr7\t116771990\tMET:c.3028+1G>T\tG\tT\t100\tPASS\t"
            "GENE=MET;RNA_EFFECT=exon_14_skipping\n"
        )


def write_structural_variant_bam(path: Path) -> None:
    """Write evidence for DEL, tandem-DUP, INV, and interchromosomal TRA."""
    rng = Random(31_041)
    header_dict = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [
            {"SN": "chr1", "LN": 10_000},
            {"SN": "chr2", "LN": 10_000},
        ],
        "RG": [{"ID": "SV", "SM": "Tumour_SV_demo"}],
    }
    header = pysam.AlignmentHeader.from_dict(header_dict)
    reads = []

    def make_read(
        name: str, start: int, cigar, flag: int = 0, mapq: int = 60,
        reference_id: int = 0,
    ):
        query_length = sum(
            length for operation, length in cigar
            if operation in (0, 1, 4, 7, 8)
        )
        read = pysam.AlignedSegment(header)
        read.query_name = name
        read.query_sequence = ("ACGT" * ((query_length + 3) // 4))[:query_length]
        read.flag = flag
        read.reference_id = reference_id
        read.reference_start = start
        read.mapping_quality = mapq
        read.cigar = cigar
        read.query_qualities = pysam.qualitystring_to_array("I" * query_length)
        read.set_tag("RG", "SV")
        return read

    def add_pair(
        name: str, left_start: int, right_start: int,
        orientation: str = "FR", proper: bool = False, mapq: int = 60,
    ) -> None:
        orientations = {
            "FR": (False, True),
            "RF": (True, False),
            "FF": (False, False),
            "RR": (True, True),
        }
        left_reverse, right_reverse = orientations[orientation]
        fragment_length = right_start + 100 - left_start
        proper_bit = 2 if proper else 0
        left_flag = 1 | proper_bit | 64
        right_flag = 1 | proper_bit | 128
        if left_reverse:
            left_flag |= 16
            right_flag |= 32
        if right_reverse:
            right_flag |= 16
            left_flag |= 32
        left = make_read(
            name, left_start, [(0, 100)], flag=left_flag, mapq=mapq
        )
        right = make_read(
            name, right_start, [(0, 100)], flag=right_flag, mapq=mapq
        )
        left.next_reference_id = 0
        right.next_reference_id = 0
        left.next_reference_start = right_start
        right.next_reference_start = left_start
        left.template_length = fragment_length
        right.template_length = -fragment_length
        reads.extend((left, right))

    def add_interchrom_pair(
        name: str, chr1_start: int, chr2_start: int,
    ) -> None:
        left = make_read(name, chr1_start, [(0, 100)], flag=1 | 64 | 32)
        right = make_read(
            name, chr2_start, [(0, 100)], flag=1 | 128 | 16,
            reference_id=1,
        )
        left.next_reference_id = 1
        right.next_reference_id = 0
        left.next_reference_start = chr2_start
        right.next_reference_start = chr1_start
        left.template_length = 0
        right.template_length = 0
        reads.extend((left, right))

    def add_softclips(
        prefix: str, left_breakpoint: int, right_breakpoint: int,
        reference_id: int = 0,
    ) -> None:
        for read_index in range(6):
            reads.append(make_read(
                f"{prefix}_left_clip_{read_index + 1:02d}", left_breakpoint - 80,
                [(0, 80), (4, 40)], flag=16 if read_index % 2 else 0,
                reference_id=reference_id,
            ))
            reads.append(make_read(
                f"{prefix}_right_clip_{read_index + 1:02d}", right_breakpoint,
                [(4, 40), (0, 80)], flag=16 if read_index % 2 else 0,
                reference_id=reference_id,
            ))

    def cigar_string(cigar) -> str:
        operation_names = {0: "M", 1: "I", 2: "D", 3: "N", 4: "S"}
        return "".join(f"{length}{operation_names[operation]}" for operation, length in cigar)

    def add_split_reads(
        prefix: str, left_start: int, left_cigar,
        right_start: int, right_cigar, right_reverse: bool = False,
        left_reference_id: int = 0, right_reference_id: int = 0,
    ) -> None:
        left_chrom = header.references[left_reference_id]
        right_chrom = header.references[right_reference_id]
        right_strand = "-" if right_reverse else "+"
        for read_index in range(6):
            name = f"{prefix}_split_{read_index + 1:02d}"
            left = make_read(
                name, left_start, left_cigar, reference_id=left_reference_id
            )
            right_flag = 2048 | (16 if right_reverse else 0)
            right = make_read(
                name, right_start, right_cigar, flag=right_flag,
                reference_id=right_reference_id,
            )
            left.set_tag(
                "SA",
                f"{right_chrom},{right_start + 1},{right_strand},"
                f"{cigar_string(right_cigar)},60,0;",
            )
            right.set_tag(
                "SA", f"{left_chrom},{left_start + 1},+,"
                f"{cigar_string(left_cigar)},60,0;"
            )
            reads.extend((left, right))

    # Dense background reads span the complete simulated locus. Inside the
    # heterozygous deletion one background component is absent, but the normal
    # fragment component below remains, producing a partial rather than empty
    # depth loss.
    background_index = 0
    while background_index < 220:
        start = rng.randint(1_050, 8_350)
        overlaps_deletion = start < 3_000 and start + 100 > 2_000
        if overlaps_deletion and rng.random() < 0.52:
            continue
        background_index += 1
        reads.append(make_read(
            f"background_{background_index:03d}", start, [(0, 100)],
            flag=16 if rng.random() < 0.5 else 0,
            mapq=rng.randint(43, 60),
        ))

    # Extra depth over the tandem duplication produces the expected gain.
    for read_index in range(1, 37):
        start = rng.randint(3_950, 4_750)
        reads.append(make_read(
            f"duplication_depth_{read_index:02d}", start, [(0, 100)],
            flag=16 if rng.random() < 0.5 else 0,
            mapq=rng.randint(48, 60),
        ))

    # Concordant FR fragments continue across every SV locus. Besides making
    # the mixed-sample background realistic, they establish a tight insert-size
    # baseline for classifying the event-supporting pairs.
    pair_index = 0
    while pair_index < 103:
        left_start = rng.randint(1_050, 8_050)
        insert_size = max(155, min(235, round(rng.gauss(195, 16))))
        right_start = left_start + insert_size
        overlaps_deletion = left_start < 3_000 and right_start + 100 > 2_000
        if overlaps_deletion and rng.random() < 0.48:
            continue
        pair_index += 1
        add_pair(
            f"normal_pair_{pair_index:03d}", left_start,
            right_start, orientation="FR", proper=True,
            mapq=rng.randint(45, 60),
        )

    # DEL: direct CIGAR gaps, large-insert FR pairs, split reads, and clips.
    for read_index in range(12):
        left_match = 52 + (read_index * 7) % 17
        right_match = 120 - left_match
        reads.append(make_read(
            f"del_cigar_{read_index + 1:02d}", 2_000 - left_match,
            [(0, left_match), (2, 1_000), (0, right_match)],
            flag=16 if read_index % 2 else 0,
        ))
    for pair_index in range(8):
        add_pair(
            f"del_large_insert_{pair_index + 1:02d}",
            1_740 + pair_index * 10, 3_070 + pair_index * 10,
            orientation="FR",
        )
    add_softclips("del", 2_000, 3_000)
    add_split_reads(
        "del", 1_930, [(0, 70), (4, 50)],
        3_000, [(4, 50), (0, 70)],
    )

    # DUP: RF/everted pairs point outwards across the tandem junction.
    for pair_index in range(8):
        add_pair(
            f"dup_everted_{pair_index + 1:02d}",
            3_820 + pair_index * 9, 4_820 + pair_index * 9,
            orientation="RF",
        )
    add_softclips("dup", 4_000, 4_800)
    add_split_reads(
        "dup", 4_730, [(0, 70), (4, 50)],
        4_000, [(4, 50), (0, 70)],
    )

    # INV: same-strand FF and RR pairs bracket the inverted segment.
    for pair_index in range(6):
        add_pair(
            f"inv_ff_{pair_index + 1:02d}",
            5_780 + pair_index * 10, 6_920 + pair_index * 10,
            orientation="FF",
        )
        add_pair(
            f"inv_rr_{pair_index + 1:02d}",
            5_850 + pair_index * 10, 6_990 + pair_index * 10,
            orientation="RR",
        )
    add_softclips("inv", 6_000, 6_900)
    add_split_reads(
        "inv", 5_930, [(0, 70), (4, 50)],
        6_900, [(4, 50), (0, 70)], right_reverse=True,
    )

    # TRA: reciprocal chr1/chr2 pairs, clips, and chimeric split alignments.
    for pair_index in range(12):
        add_interchrom_pair(
            f"tra_chr1_chr2_{pair_index + 1:02d}",
            7_470 + pair_index * 8, 4_950 + pair_index * 6,
        )
    add_softclips("tra_chr1", 7_600, 7_600)
    add_softclips("tra_chr2", 5_000, 5_000, reference_id=1)
    add_split_reads(
        "tra", 7_530, [(0, 70), (4, 50)],
        5_000, [(4, 50), (0, 70)],
        left_reference_id=0, right_reference_id=1,
    )

    # Local chr2 coverage makes the reciprocal breakpoint useful in mate views.
    for read_index, start in enumerate(range(4_300, 5_581, 32), start=1):
        reads.append(make_read(
            f"chr2_background_{read_index:02d}", start, [(0, 100)],
            flag=16 if read_index % 2 else 0, mapq=55,
            reference_id=1,
        ))

    reads.sort(key=lambda read: (
        read.reference_id, read.reference_start, read.query_name, read.flag,
    ))
    with pysam.AlignmentFile(str(path), "wb", header=header) as bam:
        for read in reads:
            bam.write(read)
    index_path = Path(f"{path}.bai")
    if index_path.exists():
        index_path.unlink()
    pysam.index(str(path))


def write_structural_variant_vcf(path: Path) -> None:
    """Write DEL, tandem-DUP, INV, and reciprocal interchromosomal BNDs."""
    with path.open("w", encoding="utf-8") as handle:
        handle.write("##fileformat=VCFv4.2\n")
        handle.write("##contig=<ID=chr1,length=10000>\n")
        handle.write("##contig=<ID=chr2,length=10000>\n")
        handle.write(
            '##INFO=<ID=END,Number=1,Type=Integer,Description="End coordinate">\n'
        )
        handle.write(
            '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV type">\n'
        )
        handle.write(
            '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="SV length">\n'
        )
        handle.write(
            '##INFO=<ID=MATEID,Number=1,Type=String,Description="ID of mate breakend">\n'
        )
        handle.write(
            '##INFO=<ID=EVENT,Number=1,Type=String,Description="Breakend event ID">\n'
        )
        handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        handle.write(
            "chr1\t2000\tDEL1000\tN\t<DEL>\t100\tPASS\t"
            "END=3000;SVTYPE=DEL;SVLEN=-1000\n"
        )
        handle.write(
            "chr1\t4000\tDUP800\tN\t<DUP:TANDEM>\t100\tPASS\t"
            "END=4800;SVTYPE=DUP;SVLEN=800\n"
        )
        handle.write(
            "chr1\t6000\tINV900\tN\t<INV>\t100\tPASS\t"
            "END=6900;SVTYPE=INV;SVLEN=900\n"
        )
        handle.write(
            "chr1\t7600\tTRA_chr1_chr2_A\tN\tN]chr2:5001]\t100\tPASS\t"
            "SVTYPE=BND;MATEID=TRA_chr1_chr2_B;EVENT=TRA1\n"
        )
        handle.write(
            "chr2\t5001\tTRA_chr1_chr2_B\tN\t[chr1:7600[N\t100\tPASS\t"
            "SVTYPE=BND;MATEID=TRA_chr1_chr2_A;EVENT=TRA1\n"
        )


def refresh_tabix(path: Path, preset: str) -> None:
    compressed_path = Path(f"{path}.gz")
    index_path = Path(f"{compressed_path}.tbi")
    if compressed_path.exists():
        compressed_path.unlink()
    if index_path.exists():
        index_path.unlink()
    pysam.tabix_compress(str(path), str(compressed_path), force=True)
    pysam.tabix_index(str(compressed_path), preset=preset, force=True)


def write_chip_signal(path: Path, baseline: float, peaks) -> None:
    window_start = 101_865_500
    window_end = 101_869_500
    bin_width = 10
    with path.open("w", encoding="utf-8") as handle:
        bin_index = 0
        position = window_start
        while position < window_end:
            bin_end = min(position + bin_width, window_end)
            midpoint = (position + bin_end) / 2
            value = baseline + ((bin_index * 17) % 11) * 0.025
            for center, amplitude, width in peaks:
                distance = midpoint - center
                value += amplitude * exp(
                    -(distance * distance) / (2 * width * width)
                )
            handle.write(f"chr9\t{position}\t{bin_end}\t{value:.3f}\n")
            position = bin_end
            bin_index += 1


def write_multi_sample_review_regions(path: Path) -> None:
    """Candidate windows used by the reproducible multi-sample HTML report."""
    path.write_text(
        "# chrom\tstart\tend\tname\n"
        "chrDemo\t78\t119\ttumour_subclonal_snv\n"
        "chrDemo\t98\t139\tshared_heterozygous_snv\n"
        "chrDemo\t139\t180\trelapse_enriched_snv\n",
        encoding="utf-8",
    )


def write_hic_tracks(tad_path: Path, loop_path: Path, contact_path: Path) -> None:
    """Write nested TAD calls, loop calls, and a binned contact map."""
    tad_path.write_text(
        "chrDemo\t0\t125\tTAD-A\t0.78\n"
        "chrDemo\t125\t275\tTAD-B\t0.96\n"
        "chrDemo\t150\t230\tSubTAD-B1\t0.64\n"
        "chrDemo\t275\t400\tTAD-C\t0.83\n",
        encoding="utf-8",
    )
    loop_path.write_text(
        "chrDemo\t18\t28\tchrDemo\t102\t112\tA-boundary\t18\n"
        "chrDemo\t62\t72\tchrDemo\t172\t182\tA-B contact\t9\n"
        "chrDemo\t132\t142\tchrDemo\t258\t268\tB-boundary\t30\n"
        "chrDemo\t158\t168\tchrDemo\t218\t228\tSubTAD loop\t24\n"
        "chrDemo\t286\t296\tchrDemo\t374\t384\tC-boundary\t16\n",
        encoding="utf-8",
    )
    loop_boosts = {
        (1, 5): 32.0, (3, 8): 18.0, (7, 13): 40.0,
        (8, 11): 34.0, (14, 19): 26.0,
    }
    contact_lines = []
    bin_size = 20
    bin_count = 20
    for first_bin in range(bin_count):
        first_midpoint = first_bin * bin_size + bin_size / 2
        first_domain = 0 if first_midpoint < 125 else (1 if first_midpoint < 275 else 2)
        for second_bin in range(first_bin, bin_count):
            second_midpoint = second_bin * bin_size + bin_size / 2
            second_domain = 0 if second_midpoint < 125 else (1 if second_midpoint < 275 else 2)
            distance = second_bin - first_bin
            score = 82.0 * exp(-distance / 2.7)
            score *= 1.30 if first_domain == second_domain else 0.42
            score += loop_boosts.get((first_bin, second_bin), 0.0)
            score += ((first_bin * 7 + second_bin * 11) % 5) * 0.7
            start1 = first_bin * bin_size
            start2 = second_bin * bin_size
            contact_lines.append(
                f"chrDemo\t{start1}\t{start1 + bin_size}\t"
                f"chrDemo\t{start2}\t{start2 + bin_size}\t.\t{score:.3f}\n"
            )
    contact_path.write_text("".join(contact_lines), encoding="utf-8")


def main() -> None:
    ensure_demo_directories()
    tumour_profile = {
        95: 0.28, 104: 0.36, 118: 0.48, 132: 0.31,
        145: 0.57, 158: 0.27, 169: 0.41,
    }
    normal_profile = {118: 0.49, 145: 0.06}
    relapse_profile = {
        104: 0.24, 118: 0.51, 145: 0.64, 158: 0.43, 169: 0.33,
    }
    small_reference = write_small_reference(REFERENCE_PATH)
    long_reference = write_long_reference(LONG_REFERENCE_PATH)
    rna_references = write_rna_fusion_reference(RNA_REFERENCE_PATH)
    write_rna_fusion_gtf(ANNOTATIONS_DIR / "demo_rna_fusion.gtf")
    write_rna_fusion_bam(
        ALIGNMENTS_DIR / "demo_rna_fusion.bam", rna_references
    )
    write_eml4_alk_dna_bam(
        ALIGNMENTS_DIR / "demo_eml4_alk_dna.bam", rna_references
    )
    write_long_read_bam(
        ALIGNMENTS_DIR / "demo_long_reads.bam", long_reference
    )
    write_small_variant_vcf(
        VARIANTS_DIR / "demo_tumour.vcf", small_reference,
        tumour_profile, "T",
    )
    write_small_variant_vcf(
        VARIANTS_DIR / "demo_relapse.vcf", small_reference,
        relapse_profile, "R",
    )
    write_variant_bam(
        ALIGNMENTS_DIR / "demo_tumour.bam", "Tumour", 240,
        tumour_profile, seed=9_601,
    )
    write_variant_bam(
        ALIGNMENTS_DIR / "demo_normal.bam", "Normal", 180,
        normal_profile, seed=7_201,
    )
    write_variant_bam(
        ALIGNMENTS_DIR / "demo_relapse.bam", "Relapse", 210,
        relapse_profile, seed=8_401,
    )
    write_variant_bam(
        ALIGNMENTS_DIR / "demo_tagged_reads.bam", "Tumour", 210,
        tumour_profile, seed=9_701,
        read_groups=("Library_A", "Library_B", "Library_C"),
    )
    write_molecule_bam(
        ALIGNMENTS_DIR / "demo_molecule_reads.bam", small_reference
    )
    write_insertion_bam(ALIGNMENTS_DIR / "demo_insertions.bam")
    cnv_reference = write_cnv_reference(CNV_REFERENCE_PATH)
    cnv_sites = build_cnv_variant_sites(cnv_reference)
    cnv_read_count = write_cnv_tumour_bam(
        ALIGNMENTS_DIR / "demo_cnv_tumour.bam", cnv_reference, cnv_sites,
    )
    write_cnv_baf_vcf(VARIANTS_DIR / "demo_baf.vcf", cnv_sites)
    write_copy_number_tracks(
        ANNOTATIONS_DIR / "demo_cnv.seg",
        ANNOTATIONS_DIR / "demo_cnv_states.bed",
    )
    write_met_ex14_bam(ALIGNMENTS_DIR / "demo_met_ex14.bam")
    write_met_ex14_tracks()
    write_structural_variant_bam(ALIGNMENTS_DIR / "demo_structural_variants.bam")
    write_structural_variant_vcf(VARIANTS_DIR / "demo_structural_variants.vcf")
    write_multi_sample_review_regions(
        ANNOTATIONS_DIR / "demo_multi_sample_review.bed"
    )
    write_hic_tracks(
        ANNOTATIONS_DIR / "demo_hic_domains.tad",
        ANNOTATIONS_DIR / "demo_hic_loops.bedpe",
        ANNOTATIONS_DIR / "demo_hic_contacts.bedpe",
    )

    write_chip_signal(
        SIGNALS_DIR / "demo_ctcf_control.signal", 0.35,
        [(101_866_220, 48.0, 65), (101_867_360, 7.0, 110),
         (101_868_650, 42.0, 85)],
    )
    write_chip_signal(
        SIGNALS_DIR / "demo_ctcf_knockdown.signal", 0.25,
        [(101_866_220, 13.0, 70), (101_867_360, 3.5, 120),
         (101_868_650, 10.0, 90)],
    )
    write_chip_signal(
        SIGNALS_DIR / "demo_ctcf_mel.signal", 0.30,
        [(101_866_220, 55.0, 115), (101_867_250, 18.0, 190),
         (101_868_650, 30.0, 135)],
    )

    refresh_tabix(VARIANTS_DIR / "demo_variants.vcf", "vcf")
    refresh_tabix(VARIANTS_DIR / "demo_tumour.vcf", "vcf")
    refresh_tabix(VARIANTS_DIR / "demo_relapse.vcf", "vcf")
    refresh_tabix(VARIANTS_DIR / "demo_baf.vcf", "vcf")
    refresh_tabix(VARIANTS_DIR / "demo_met_ex14.vcf", "vcf")
    refresh_tabix(VARIANTS_DIR / "demo_structural_variants.vcf", "vcf")
    refresh_tabix(ANNOTATIONS_DIR / "demo_dnase.narrowPeak", "bed")
    refresh_tabix(SIGNALS_DIR / "demo_ctcf_control.signal", "bed")
    refresh_tabix(SIGNALS_DIR / "demo_ctcf_knockdown.signal", "bed")
    refresh_tabix(SIGNALS_DIR / "demo_ctcf_mel.signal", "bed")
    print(
        "Regenerated demo inputs in out/demo_data/ "
        f"({cnv_read_count:,} CNV reads; {len(cnv_sites)} BAF loci)."
    )


if __name__ == "__main__":
    main()
