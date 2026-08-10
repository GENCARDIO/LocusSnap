import os
import sys

import pysam

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from locus_snap.read_model import fetch_reads, open_alignment_file
from locus_snap.snapshot import BamSnapshot

CHROM = "chr1"
REF_LENGTH = 2000
N_READS = 5


def _build_bam_and_cram(tmp_path):
    sequence = ("ACGT" * (REF_LENGTH // 4 + 1))[:REF_LENGTH]

    fasta_path = tmp_path / "reference.fa"
    with fasta_path.open("w", encoding="utf-8") as handle:
        handle.write(f">{CHROM}\n")
        for i in range(0, REF_LENGTH, 60):
            handle.write(sequence[i:i + 60] + "\n")
    pysam.faidx(str(fasta_path))

    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": CHROM, "LN": REF_LENGTH}]}
    bam_path = tmp_path / "reads.bam"
    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as bam:
        for i in range(N_READS):
            start = 100 + i * 50
            read_sequence = sequence[start:start + 50]
            segment = pysam.AlignedSegment(bam.header)
            segment.query_name = f"read{i}"
            segment.query_sequence = read_sequence
            segment.flag = 0
            segment.reference_id = 0
            segment.reference_start = start
            segment.mapping_quality = 60
            segment.cigar = [(0, len(read_sequence))]
            segment.query_qualities = pysam.qualitystring_to_array("I" * len(read_sequence))
            bam.write(segment)
    pysam.index(str(bam_path))

    cram_path = tmp_path / "reads.cram"
    with pysam.AlignmentFile(str(bam_path), "rb") as src, pysam.AlignmentFile(
        str(cram_path), "wc", template=src, reference_filename=str(fasta_path)
    ) as dst:
        for segment in src:
            dst.write(segment)
    pysam.index(str(cram_path))

    return str(bam_path), str(cram_path), str(fasta_path)


def test_open_alignment_file_selects_mode_by_extension(monkeypatch):
    calls = []

    def fake_open(path, mode, **kwargs):
        calls.append((path, mode, kwargs))
        return "handle"

    monkeypatch.setattr(pysam, "AlignmentFile", fake_open)

    assert open_alignment_file("sample.cram", reference="ref.fa") == "handle"
    assert calls[-1] == ("sample.cram", "rc", {"reference_filename": "ref.fa"})

    assert open_alignment_file("sample.bam") == "handle"
    assert calls[-1] == ("sample.bam", "rb", {})

    assert open_alignment_file("sample.CRAM", reference="ref.fa") == "handle"
    assert calls[-1] == ("sample.CRAM", "rc", {"reference_filename": "ref.fa"})


def test_fetch_reads_from_cram_matches_bam(tmp_path):
    bam_path, cram_path, fasta_path = _build_bam_and_cram(tmp_path)

    bam_reads = fetch_reads(bam_path, CHROM, 0, REF_LENGTH)
    cram_reads = fetch_reads(cram_path, CHROM, 0, REF_LENGTH)

    assert len(bam_reads) == N_READS
    assert len(cram_reads) == N_READS
    assert sorted(r.query_name for r in bam_reads) == sorted(r.query_name for r in cram_reads)


def test_bam_snapshot_renders_from_cram(tmp_path):
    _, cram_path, fasta_path = _build_bam_and_cram(tmp_path)

    snap = BamSnapshot(
        bam=cram_path, chrom=CHROM, start=0, end=REF_LENGTH, fasta=fasta_path,
        output_dir=str(tmp_path), output_name="cram_snapshot",
        show_ideogram=False, annotation_sources=[],
    )
    summary = snap.snap()

    assert summary.n_reads == N_READS
    assert os.path.isfile(snap.output_path)
