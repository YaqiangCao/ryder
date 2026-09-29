import subprocess

import pytest

from src import paw


def test_bedgraph_removed_after_successful_bigwig_conversion(tmp_path, monkeypatch):
    bedgraph = tmp_path / "normalized.bdg"
    chrom_sizes = tmp_path / "chrom.sizes"
    bigwig = tmp_path / "normalized.bw"
    bedgraph.write_text("chr1\t0\t10\t1.0\n")
    chrom_sizes.write_text("chr1\t100\n")

    def successful_conversion(command, check):
        assert command == [
            "bedGraphToBigWig",
            str(bedgraph),
            str(chrom_sizes),
            str(bigwig),
        ]
        assert check is True
        bigwig.write_bytes(b"bigWig output")

    monkeypatch.setattr(paw.subprocess, "run", successful_conversion)
    paw._convertBedGraphToBigWig(str(bedgraph), str(chrom_sizes), str(bigwig))

    assert bigwig.is_file()
    assert not bedgraph.exists()


def test_bedgraph_retained_when_conversion_fails(tmp_path, monkeypatch):
    bedgraph = tmp_path / "normalized.bdg"
    chrom_sizes = tmp_path / "chrom.sizes"
    bigwig = tmp_path / "normalized.bw"
    bedgraph.write_text("chr1\t0\t10\t1.0\n")
    chrom_sizes.write_text("chr1\t100\n")

    def failed_conversion(command, check):
        raise subprocess.CalledProcessError(1, command)

    monkeypatch.setattr(paw.subprocess, "run", failed_conversion)
    with pytest.raises(RuntimeError, match="retaining intermediate file"):
        paw._convertBedGraphToBigWig(
            str(bedgraph), str(chrom_sizes), str(bigwig)
        )

    assert bedgraph.is_file()
    assert not bigwig.exists()


def test_bedgraph_retained_when_bigwig_is_not_created(tmp_path, monkeypatch):
    bedgraph = tmp_path / "normalized.bdg"
    chrom_sizes = tmp_path / "chrom.sizes"
    bigwig = tmp_path / "normalized.bw"
    bedgraph.write_text("chr1\t0\t10\t1.0\n")
    chrom_sizes.write_text("chr1\t100\n")

    monkeypatch.setattr(paw.subprocess, "run", lambda command, check: None)
    with pytest.raises(RuntimeError, match="nonempty output"):
        paw._convertBedGraphToBigWig(
            str(bedgraph), str(chrom_sizes), str(bigwig)
        )

    assert bedgraph.is_file()
    assert not bigwig.exists()
