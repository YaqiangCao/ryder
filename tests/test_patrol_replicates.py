import gzip

import numpy as np
import pandas as pd
import pyBigWig

from src.patrol_replicates import (
    benjamini_hochberg,
    build_normalization_factors,
    count_bedpe,
    estimate_dispersions,
    fit_negative_binomial,
    read_regions,
)


def write_bigwig(path, values):
    with pyBigWig.open(str(path), "w") as bw:
        bw.addHeader([("chr1", 1_000)])
        starts = [100, 300]
        bw.addEntries(["chr1", "chr1"], starts, ends=[200, 400], values=values)


def test_bedpe_fragment_counting(tmp_path):
    bed = tmp_path / "regions.bed"
    bed.write_text("chr1\t100\t200\nchr1\t300\t400\n")
    bedpe = tmp_path / "reads.bedpe.gz"
    rows = (
        "chr1\t90\t120\tchr1\t150\t180\ta\t30\t+\t-\n"
        "chr1\t250\t280\tchr1\t320\t350\tb\t30\t+\t-\n"
        "chr1\t120\t150\tchr1\t330\t360\tc\t30\t+\t-\n"
        "chr1\t120\t150\tchr1\t330\t360\td\t5\t+\t-\n"
    )
    with gzip.open(bedpe, "wt") as handle:
        handle.write(rows)
    counts, library_size, malformed = count_bedpe(
        str(bedpe), read_regions(str(bed)), mapq_cutoff=10
    )
    np.testing.assert_array_equal(counts, [2, 2])
    assert library_size == 3
    assert malformed == 0


def test_paw_multiplier_is_converted_to_inverse_exposure(tmp_path):
    bed = tmp_path / "regions.bed"
    bed.write_text("chr1\t100\t200\nchr1\t300\t400\n")
    regions = read_regions(str(bed))
    raw1, raw2 = tmp_path / "raw1.bw", tmp_path / "raw2.bw"
    norm1, norm2 = tmp_path / "norm1.bw", tmp_path / "norm2.bw"
    write_bigwig(raw1, [1.0, 2.0])
    write_bigwig(norm1, [1.0, 2.0])
    write_bigwig(raw2, [1.0, 2.0])
    write_bigwig(norm2, [2.0, 4.0])
    samples = pd.DataFrame(
        {
            "sample": ["a", "b"],
            "raw_bw": [str(raw1), str(raw2)],
            "normalized_bw": [str(norm1), str(norm2)],
        }
    )
    factors, multipliers, _ = build_normalization_factors(
        regions, samples, np.array([100.0, 200.0])
    )
    np.testing.assert_allclose(multipliers, [[1, 2], [1, 2]])
    np.testing.assert_allclose(factors, np.ones((2, 2)))


def test_bh_and_negative_binomial_effect_direction():
    stable = np.tile(np.array([100, 105, 98, 102]), (40, 1))
    stable += np.arange(40)[:, None] % 7
    changed = np.array([[100, 110, 800, 850]])
    counts = np.vstack([stable, changed]).astype(float)
    factors = np.ones_like(counts)
    conditions = np.array(["WT", "WT", "KO", "KO"])
    _, _, dispersions, _ = estimate_dispersions(counts, factors, conditions)
    log2fc, _, _, pvalue, padj = fit_negative_binomial(
        counts, factors, conditions, "WT", "KO", dispersions, min_count=10
    )
    assert log2fc[-1] > 2
    assert pvalue[-1] < 0.01
    np.testing.assert_allclose(padj, benjamini_hochberg(pvalue), equal_nan=True)
