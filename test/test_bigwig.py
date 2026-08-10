import os
import sys

import pyBigWig
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from locus_snap.annotations import AnnotationSource, infer_track_format
from locus_snap.cli import main

TEST_BAM = os.path.join(os.path.dirname(__file__), "test.bam")


def write_bigwig(path, chrom="chr1", length=1000, entries=()):
    bw = pyBigWig.open(str(path), "w")
    bw.addHeader([(chrom, length)])
    if entries:
        starts = [entry[0] for entry in entries]
        ends = [entry[1] for entry in entries]
        values = [entry[2] for entry in entries]
        bw.addEntries([chrom] * len(entries), starts, ends=ends, values=values)
    bw.close()


def test_format_inference_maps_bigwig_extensions_to_signal():
    assert infer_track_format("sample.bw") == "signal"
    assert infer_track_format("sample.bigWig") == "signal"
    assert infer_track_format("sample.bw.gz") == "signal"


def test_bigwig_track_loads_as_signal_kind(tmp_path):
    bw_path = tmp_path / "coverage.bw"
    write_bigwig(bw_path, entries=[(0, 50, 1.5), (100, 200, 2.5)])

    source = AnnotationSource(str(bw_path))
    assert source.kind == "signal"
    assert source.is_bigwig is True

    track = source.fetch("chr1", 0, 250)
    assert track.kind == "signal"
    assert [item.value for item in track.items] == [1.5, 2.5]
    assert [(item.start, item.end) for item in track.items] == [(0, 50), (100, 200)]
    assert track.rows == [track.items]


def test_bigwig_fetch_clips_to_requested_window(tmp_path):
    bw_path = tmp_path / "coverage.bw"
    write_bigwig(bw_path, entries=[(0, 300, 4.0)])

    track = AnnotationSource(str(bw_path)).fetch("chr1", 100, 200)
    assert [(item.start, item.end) for item in track.items] == [(100, 200)]


def test_bigwig_fetch_resolves_chr_prefix_mismatch(tmp_path):
    bw_path = tmp_path / "coverage.bw"
    write_bigwig(bw_path, chrom="1", entries=[(0, 50, 1.0)])

    track = AnnotationSource(str(bw_path)).fetch("chr1", 0, 100)
    assert [item.value for item in track.items] == [1.0]


def test_bigwig_fetch_returns_empty_for_unknown_contig(tmp_path):
    bw_path = tmp_path / "coverage.bw"
    write_bigwig(bw_path, chrom="chr1", entries=[(0, 50, 1.0)])

    track = AnnotationSource(str(bw_path)).fetch("chr9", 0, 100)
    assert track.items == []


def test_bigwig_signal_kind_rejects_negative_values(tmp_path):
    bw_path = tmp_path / "signed.bw"
    write_bigwig(bw_path, entries=[(0, 50, -1.0)])

    with pytest.raises(ValueError, match="require non-negative"):
        AnnotationSource(str(bw_path)).fetch("chr1", 0, 100)


def test_bigwig_explicit_log2_kind_allows_negative_values(tmp_path):
    bw_path = tmp_path / "signed.bw"
    write_bigwig(bw_path, entries=[(0, 50, -1.0), (100, 150, 2.0)])

    source = AnnotationSource(str(bw_path), kind="log2")
    track = source.fetch("chr1", 0, 200)
    assert [item.value for item in track.items] == [-1.0, 2.0]


def test_bigwig_type_alias_accepted_in_custom_track_kind(tmp_path):
    bw_path = tmp_path / "coverage.bw"
    write_bigwig(bw_path, entries=[(0, 50, 1.0)])

    assert AnnotationSource(str(bw_path), kind="bigwig").kind == "signal"
    assert AnnotationSource(str(bw_path), kind="bw").kind == "signal"


def test_invalid_bigwig_file_is_rejected_at_construction(tmp_path):
    bad_path = tmp_path / "not_really.bw"
    bad_path.write_text("this is not a bigwig file", encoding="utf-8")

    with pytest.raises(ValueError, match="Cannot open BigWig"):
        AnnotationSource(str(bad_path))


def test_bigwig_does_not_require_tabix_index(tmp_path):
    bw_path = tmp_path / "coverage.bw"
    write_bigwig(bw_path, entries=[(0, 50, 1.0)])

    # No .tbi/.csi alongside it - must not raise the compressed-track index error.
    source = AnnotationSource(str(bw_path))
    assert source.compressed is False


def test_cli_renders_with_a_bigwig_track(tmp_path):
    bw_path = tmp_path / "coverage.bw"
    write_bigwig(
        bw_path, chrom="chr9", length=141213431,
        entries=[(101867470, 101867630, 5.0)],
    )

    rc = main([
        "--bam", TEST_BAM,
        "--region", "chr9:101867481-101867620",
        "--track", str(bw_path),
        "--refseq", "none",
        "--output_dir", str(tmp_path),
        "--output_name", "with_bigwig",
    ])

    assert rc == 0
    assert (tmp_path / "with_bigwig.png").is_file()
