#!/usr/bin/env python3
# --coding:utf-8--
"""Replicate-aware differential-region testing after PAW normalization.

The command keeps raw fragment counts as the response variable and represents
PAW's potentially region-specific correction as a normalization-factor matrix
in a negative-binomial GLM.  This follows the DESeq2 model

    K_ij ~ NB(mu_ij, alpha_i),  mu_ij = NF_ij q_ij

without rounding continuous normalized bigWig signals into pseudo-counts.
It is a DESeq2-style implementation, not a bit-for-bit reimplementation of
the DESeq2 or PyDESeq2 dispersion estimators.
"""

from __future__ import annotations

import gzip
import json
import math
import os
from collections import defaultdict
from datetime import datetime
from pathlib import Path

import click
import matplotlib as mpl

mpl.use("pdf")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyBigWig
import statsmodels.api as sm
from scipy.optimize import least_squares
from scipy.stats import norm

mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["savefig.bbox"] = "tight"
mpl.rcParams["savefig.transparent"] = True
mpl.rcParams["font.size"] = 8.0


def rprint(message: str) -> None:
    print(f"{datetime.now()}\t{message}", flush=True)


def read_regions(filepath: str) -> pd.DataFrame:
    """Read unique BED intervals of at least 100 bp, matching PATROL policy."""
    records = []
    seen = set()
    with open(filepath) as handle:
        for line_number, line in enumerate(handle, 1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            try:
                start, end = int(fields[1]), int(fields[2])
            except ValueError:
                continue
            if start < 0 or end - start < 100:
                continue
            region_id = f"{fields[0]}:{start}-{end}"
            if region_id in seen:
                continue
            seen.add(region_id)
            records.append((region_id, fields[0], start, end, line_number))
    if not records:
        raise ValueError(f"No valid regions of at least 100 bp in {filepath}")
    return pd.DataFrame(
        records, columns=["region", "chrom", "start", "end", "bed_line"]
    ).set_index("region")


def read_sample_sheet(filepath: str) -> pd.DataFrame:
    """Read and validate sample, condition, BEDPE, raw-bw, normalized-bw data."""
    samples = pd.read_csv(filepath, sep="\t", dtype=str)
    required = ["sample", "condition", "bedpe", "raw_bw", "normalized_bw"]
    missing = [column for column in required if column not in samples.columns]
    if missing:
        raise ValueError(f"Sample sheet is missing columns: {', '.join(missing)}")
    if samples.empty or samples["sample"].duplicated().any():
        raise ValueError("Sample names must be present and unique.")
    base = Path(filepath).resolve().parent
    for column in ("bedpe", "raw_bw", "normalized_bw"):
        samples[column] = samples[column].map(
            lambda value: str((base / value).resolve())
            if not Path(value).is_absolute()
            else str(Path(value).resolve())
        )
        absent = [value for value in samples[column] if not Path(value).is_file()]
        if absent:
            raise FileNotFoundError(f"Missing {column} file(s): {absent}")
    return samples


def build_bin_index(regions: pd.DataFrame, bin_size: int = 100_000):
    """Build a small genomic-bin lookup for streaming BEDPE overlap counts."""
    index = defaultdict(list)
    starts = regions["start"].to_numpy(dtype=np.int64)
    ends = regions["end"].to_numpy(dtype=np.int64)
    chroms = regions["chrom"].to_numpy()
    for region_index, (chrom, start, end) in enumerate(zip(chroms, starts, ends)):
        for genomic_bin in range(start // bin_size, (end - 1) // bin_size + 1):
            index[(chrom, genomic_bin)].append(region_index)
    return index, starts, ends


def _overlap_ids(chrom, start, end, bin_index, starts, ends, bin_size):
    candidates = set()
    for genomic_bin in range(start // bin_size, (end - 1) // bin_size + 1):
        candidates.update(bin_index.get((chrom, genomic_bin), ()))
    return {idx for idx in candidates if starts[idx] < end and ends[idx] > start}


def count_bedpe(
    filepath: str,
    regions: pd.DataFrame,
    mapq_cutoff: float = 10,
    bin_size: int = 100_000,
):
    """Count passing paired-end fragments overlapping each candidate region."""
    bin_index, starts, ends = build_bin_index(regions, bin_size=bin_size)
    counts = np.zeros(len(regions), dtype=np.int64)
    passing_fragments = 0
    malformed = 0
    opener = gzip.open if str(filepath).endswith(".gz") else open
    with opener(filepath, "rt") as handle:
        for line_number, line in enumerate(handle, 1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                malformed += 1
                continue
            try:
                chrom1, start1, end1 = fields[0], int(fields[1]), int(fields[2])
                chrom2, start2, end2 = fields[3], int(fields[4]), int(fields[5])
                mapq = float(fields[7])
            except ValueError:
                malformed += 1
                continue
            if mapq < mapq_cutoff:
                continue
            passing_fragments += 1
            if chrom1 == chrom2:
                hit_ids = _overlap_ids(
                    chrom1,
                    min(start1, start2),
                    max(end1, end2),
                    bin_index,
                    starts,
                    ends,
                    bin_size,
                )
            else:
                hit_ids = _overlap_ids(
                    chrom1, start1, end1, bin_index, starts, ends, bin_size
                )
                hit_ids.update(
                    _overlap_ids(
                        chrom2, start2, end2, bin_index, starts, ends, bin_size
                    )
                )
            for region_index in hit_ids:
                counts[region_index] += 1
            if line_number % 5_000_000 == 0:
                rprint(f"[{Path(filepath).name}] processed {line_number:,} BEDPE rows")
    return counts, passing_fragments, malformed


def quantify_bigwig(regions: pd.DataFrame, filepath: str) -> np.ndarray:
    """Get exact integrated bigWig signal for each interval."""
    values = np.zeros(len(regions), dtype=float)
    with pyBigWig.open(filepath) as bw:
        chrom_sizes = bw.chroms()
        for index, row in enumerate(regions.itertuples()):
            if row.chrom not in chrom_sizes or row.start >= chrom_sizes[row.chrom]:
                continue
            end = min(row.end, chrom_sizes[row.chrom])
            try:
                value = bw.stats(row.chrom, row.start, end, type="sum", exact=True)[0]
            except RuntimeError:
                value = None
            if value is not None and np.isfinite(value):
                values[index] = max(float(value), 0.0)
    return values


def build_normalization_factors(
    regions: pd.DataFrame,
    samples: pd.DataFrame,
    library_sizes: np.ndarray,
):
    """Combine library depth and regional PAW correction into DESeq-style NF."""
    multipliers = np.ones((len(regions), len(samples)), dtype=float)
    diagnostics = []
    for sample_index, sample in samples.iterrows():
        raw = quantify_bigwig(regions, sample["raw_bw"])
        normalized = quantify_bigwig(regions, sample["normalized_bw"])
        valid = (raw > 0) & (normalized > 0)
        if not np.any(valid):
            raise ValueError(f"No positive raw/normalized signals for {sample['sample']}")
        ratio = np.divide(normalized, raw, out=np.full_like(raw, np.nan), where=valid)
        fallback = float(np.nanmedian(ratio))
        ratio[~np.isfinite(ratio) | (ratio <= 0)] = fallback
        multipliers[:, sample_index] = ratio
        diagnostics.append(
            {
                "sample": sample["sample"],
                "library_size": int(library_sizes[sample_index]),
                "positive_bigwig_regions": int(valid.sum()),
                "median_paw_multiplier": fallback,
            }
        )
    factors = library_sizes[np.newaxis, :] / multipliers
    row_geomean = np.exp(np.mean(np.log(factors), axis=1))
    factors /= row_geomean[:, np.newaxis]
    return factors, multipliers, pd.DataFrame(diagnostics)


def estimate_dispersions(counts: np.ndarray, factors: np.ndarray, groups: np.ndarray):
    """Estimate gene-wise dispersions, a DESeq2-form trend, and MAP-like shrinkage."""
    n_samples = counts.shape[1]
    unique_groups = np.unique(groups)
    fitted = np.zeros_like(counts, dtype=float)
    for group in unique_groups:
        mask = groups == group
        abundance = counts[:, mask].sum(axis=1) / factors[:, mask].sum(axis=1)
        fitted[:, mask] = factors[:, mask] * abundance[:, np.newaxis]
    numerator = np.sum((counts - fitted) ** 2 - fitted, axis=1)
    denominator = np.sum(fitted**2, axis=1)
    residual_df = max(n_samples - len(unique_groups), 1)
    raw_unclipped = (
        np.divide(numerator, denominator, out=np.zeros_like(numerator), where=denominator > 0)
        * n_samples
        / residual_df
    )
    raw = np.maximum(raw_unclipped, 1e-8)
    base_mean = np.mean(counts / factors, axis=1)
    # Negative method-of-moments estimates mean that observed residual variance
    # did not exceed Poisson variance. They are not valid points for fitting the
    # dispersion trend and must not be converted to a cloud at the lower bound.
    valid = (base_mean > 0) & (raw_unclipped > 0) & np.isfinite(raw_unclipped)
    if valid.sum() < 20:
        constant = float(np.median(raw[valid])) if valid.any() else 0.1
        trend = np.full_like(raw, max(constant, 1e-8))
    else:
        x = base_mean[valid]
        y = np.log(raw[valid])

        def residuals(parameters):
            predicted = np.exp(parameters[0]) + np.exp(parameters[1]) / x
            return np.log(predicted) - y

        initial = np.log([max(np.median(raw[valid]), 1e-4), 1.0])
        fit = least_squares(residuals, initial, loss="soft_l1")
        trend = np.exp(fit.x[0]) + np.exp(fit.x[1]) / np.maximum(base_mean, 1e-8)
    sampling_variance = 2.0 / residual_df
    if valid.sum() >= 2:
        log_residual = np.log(raw[valid]) - np.log(trend[valid])
        mad = np.median(np.abs(log_residual - np.median(log_residual)))
        observed_variance = (1.4826 * mad) ** 2
        prior_variance = max(observed_variance - sampling_variance, 0.25)
    else:
        # With too few positive gene-wise estimates, retain a finite weak prior
        # instead of propagating an empty-slice NaN into every GLM dispersion.
        prior_variance = 0.25
    weight = prior_variance / (prior_variance + sampling_variance)
    gene_estimate = np.where(raw_unclipped > 0, raw, trend)
    shrunk = np.exp(weight * np.log(gene_estimate) + (1 - weight) * np.log(trend))
    return raw, trend, shrunk, base_mean


def benjamini_hochberg(p_values: np.ndarray) -> np.ndarray:
    adjusted = np.full(len(p_values), np.nan)
    valid = np.isfinite(p_values)
    p = p_values[valid]
    if not len(p):
        return adjusted
    order = np.argsort(p)
    ranked = p[order]
    corrected = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    corrected = np.minimum.accumulate(corrected[::-1])[::-1]
    restored = np.empty_like(corrected)
    restored[order] = np.minimum(corrected, 1.0)
    adjusted[valid] = restored
    return adjusted


def fit_negative_binomial(
    counts: np.ndarray,
    factors: np.ndarray,
    conditions: np.ndarray,
    reference_condition: str,
    treatment_condition: str,
    dispersions: np.ndarray,
    min_count: int,
):
    """Fit one offset negative-binomial Wald test per region."""
    group = (conditions == treatment_condition).astype(float)
    design = sm.add_constant(group, prepend=True)
    log2fc = np.full(counts.shape[0], np.nan)
    standard_error = np.full(counts.shape[0], np.nan)
    statistic = np.full(counts.shape[0], np.nan)
    p_value = np.full(counts.shape[0], np.nan)
    for index, row in enumerate(counts):
        if row.sum() < min_count or np.count_nonzero(row) < 2:
            continue
        try:
            fit = sm.GLM(
                row,
                design,
                family=sm.families.NegativeBinomial(alpha=dispersions[index]),
                offset=np.log(factors[index]),
            ).fit(maxiter=100, disp=0)
        except (ValueError, np.linalg.LinAlgError, FloatingPointError):
            continue
        beta = float(fit.params[1])
        se = float(fit.bse[1])
        if not np.isfinite(beta) or not np.isfinite(se) or se <= 0:
            continue
        log2fc[index] = beta / math.log(2)
        standard_error[index] = se / math.log(2)
        statistic[index] = beta / se
        p_value[index] = 2 * norm.sf(abs(statistic[index]))
    return log2fc, standard_error, statistic, p_value, benjamini_hochberg(p_value)


def plot_results(results, counts, factors, samples, output_prefix, fdr, lfc):
    significant = (results["padj"] <= fdr) & (results["log2FoldChange"].abs() >= lfc)
    up = significant & (results["log2FoldChange"] > 0)
    down = significant & (results["log2FoldChange"] < 0)
    fig, ax = plt.subplots(figsize=(4.2, 3.2))
    ax.scatter(np.log10(results["baseMean"] + 1), results["log2FoldChange"], s=3, c="0.7")
    ax.scatter(np.log10(results.loc[down, "baseMean"] + 1), results.loc[down, "log2FoldChange"], s=5, c="#d62728", label=f"reference enriched ({down.sum():,})")
    ax.scatter(np.log10(results.loc[up, "baseMean"] + 1), results.loc[up, "log2FoldChange"], s=5, c="#2166ac", label=f"treatment enriched ({up.sum():,})")
    ax.axhline(0, color="black", linewidth=0.6)
    ax.axhline(lfc, color="0.5", linestyle="--", linewidth=0.6)
    ax.axhline(-lfc, color="0.5", linestyle="--", linewidth=0.6)
    ax.set(xlabel="log10(base mean + 1)", ylabel="log2 fold change", title=f"Replicate-aware differential regions\nFDR <= {fdr:g}")
    ax.legend(frameon=False, fontsize=7)
    fig.tight_layout()
    fig.savefig(f"{output_prefix}_MA.pdf")
    plt.close(fig)

    transformed = np.log2(counts / factors + 1).T
    transformed -= transformed.mean(axis=0)
    u, singular_values, _ = np.linalg.svd(transformed, full_matrices=False)
    coordinates = u[:, :2] * singular_values[:2]
    explained = singular_values**2 / np.sum(singular_values**2)
    fig, ax = plt.subplots(figsize=(4, 3.2))
    for condition, group_samples in samples.groupby("condition", sort=False):
        idx = group_samples.index.to_numpy()
        ax.scatter(coordinates[idx, 0], coordinates[idx, 1], s=28, label=condition)
        for position in idx:
            ax.annotate(samples.loc[position, "sample"], coordinates[position], xytext=(4, 3), textcoords="offset points", fontsize=7)
    ax.set_xlabel(f"PC1 ({explained[0] * 100:.1f}%)")
    ax.set_ylabel(f"PC2 ({explained[1] * 100:.1f}%)")
    ax.set_title("PAW-normalized region counts")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(f"{output_prefix}_PCA.pdf")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(4, 3.2))
    valid = results["baseMean"] > 0
    ax.scatter(results.loc[valid, "baseMean"], results.loc[valid, "dispersion_raw"], s=3, c="0.75", label="raw")
    order = np.argsort(results.loc[valid, "baseMean"].to_numpy())
    x = results.loc[valid, "baseMean"].to_numpy()[order]
    ax.plot(x, results.loc[valid, "dispersion_trend"].to_numpy()[order], c="#d62728", label="trend")
    ax.scatter(results.loc[valid, "baseMean"], results.loc[valid, "dispersion"], s=3, c="#2166ac", label="shrunk")
    ax.set(xscale="log", yscale="log", xlabel="base mean", ylabel="dispersion", title="Dispersion estimates")
    ax.legend(frameon=False, fontsize=7)
    fig.tight_layout()
    fig.savefig(f"{output_prefix}_dispersion.pdf")
    plt.close(fig)


def write_bed(results, mask, filepath):
    columns = ["chrom", "start", "end", "log2FoldChange", "pvalue", "padj", "baseMean"]
    results.loc[mask, columns].to_csv(filepath, sep="\t", header=False, index=False)


@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.option("-r", "regions_file", required=True, type=click.Path(exists=True), help="Candidate regions in BED format.")
@click.option("-s", "sample_sheet", required=True, type=click.Path(exists=True), help="TSV with sample, condition, bedpe, raw_bw, normalized_bw columns.")
@click.option("-o", "output_prefix", required=True, type=click.Path(), help="Output prefix.")
@click.option("--reference-condition", required=True, help="Reference condition in the sample sheet.")
@click.option("--treatment-condition", required=True, help="Treatment condition in the sample sheet.")
@click.option("--mapq", default=10.0, show_default=True, type=float, help="Minimum BEDPE MAPQ.")
@click.option("--min-count", default=10, show_default=True, type=int, help="Minimum total count required for testing.")
@click.option("--fdr", default=0.05, show_default=True, type=click.FloatRange(0, 1, min_open=True), help="BH-adjusted significance cutoff.")
@click.option("--lfc", default=1.0, show_default=True, type=float, help="Absolute log2-fold-change cutoff for BED calls.")
@click.option("--reuse-counts/--no-reuse-counts", default=False, show_default=True, help="Reuse a compatible existing count table and sample diagnostics.")
def patrol_replicates(regions_file, sample_sheet, output_prefix, reference_condition, treatment_condition, mapq, min_count, fdr, lfc, reuse_counts):
    """Test differential genomic regions using raw replicate counts and PAW offsets."""
    start_time = datetime.now()
    samples = read_sample_sheet(sample_sheet).reset_index(drop=True)
    conditions = samples["condition"].to_numpy()
    selected = {reference_condition, treatment_condition}
    if set(conditions) != selected:
        raise click.ClickException(f"Expected exactly conditions {sorted(selected)}; observed {sorted(set(conditions))}")
    condition_counts = samples["condition"].value_counts()
    if (condition_counts < 2).any():
        raise click.ClickException("At least two biological replicates per condition are required.")
    output_parent = os.path.dirname(output_prefix)
    if output_parent:
        os.makedirs(output_parent, exist_ok=True)

    regions = read_regions(regions_file)
    count_path = f"{output_prefix}_counts.tsv"
    diagnostic_path = f"{output_prefix}_sample_diagnostics.tsv"
    cached = reuse_counts and os.path.isfile(count_path) and os.path.isfile(diagnostic_path)
    if cached:
        count_table = pd.read_csv(count_path, sep="\t", index_col=0)
        previous_diagnostics = pd.read_csv(diagnostic_path, sep="\t")
        compatible = (
            count_table.index.equals(regions.index)
            and list(count_table.columns) == list(samples["sample"])
            and list(previous_diagnostics["sample"]) == list(samples["sample"])
        )
        if compatible:
            rprint(f"[{output_prefix}] reusing compatible raw count table")
            count_matrix = count_table.to_numpy(dtype=np.int64)
            library_sizes = previous_diagnostics["passing_bedpe_fragments"].to_numpy(dtype=float)
            count_diagnostics = previous_diagnostics[
                ["sample", "passing_bedpe_fragments", "malformed_bedpe_rows"]
            ].to_dict("records")
        else:
            cached = False
    if not cached:
        count_columns = []
        library_sizes = []
        count_diagnostics = []
        for sample in samples.itertuples():
            rprint(f"[{output_prefix}] counting {sample.sample}: {sample.bedpe}")
            counts, library_size, malformed = count_bedpe(sample.bedpe, regions, mapq_cutoff=mapq)
            count_columns.append(counts)
            library_sizes.append(library_size)
            count_diagnostics.append({"sample": sample.sample, "passing_bedpe_fragments": library_size, "malformed_bedpe_rows": malformed})
        count_matrix = np.column_stack(count_columns)
        library_sizes = np.asarray(library_sizes, dtype=float)
        count_table = pd.DataFrame(count_matrix, index=regions.index, columns=samples["sample"])
        count_table.to_csv(count_path, sep="\t")

    rprint(f"[{output_prefix}] deriving PAW normalization-factor matrix")
    factors, multipliers, bw_diagnostics = build_normalization_factors(regions, samples, library_sizes)
    pd.DataFrame(factors, index=regions.index, columns=samples["sample"]).to_csv(f"{output_prefix}_normalization_factors.tsv", sep="\t")
    pd.DataFrame(multipliers, index=regions.index, columns=samples["sample"]).to_csv(f"{output_prefix}_paw_multipliers.tsv", sep="\t")

    raw_dispersion, trend_dispersion, dispersion, base_mean = estimate_dispersions(count_matrix, factors, conditions)
    log2fc, lfc_se, statistic, pvalue, padj = fit_negative_binomial(
        count_matrix, factors, conditions, reference_condition, treatment_condition, dispersion, min_count
    )
    normalized_counts = count_matrix / factors
    reference_mask = conditions == reference_condition
    treatment_mask = conditions == treatment_condition
    results = regions.copy()
    results["baseMean"] = base_mean
    results[f"mean_{reference_condition}"] = normalized_counts[:, reference_mask].mean(axis=1)
    results[f"mean_{treatment_condition}"] = normalized_counts[:, treatment_mask].mean(axis=1)
    results["log2FoldChange"] = log2fc
    results["lfcSE"] = lfc_se
    results["stat"] = statistic
    results["pvalue"] = pvalue
    results["padj"] = padj
    results["dispersion_raw"] = raw_dispersion
    results["dispersion_trend"] = trend_dispersion
    results["dispersion"] = dispersion
    results.to_csv(f"{output_prefix}_results.tsv", sep="\t", index_label="region")

    significant = (results["padj"] <= fdr) & (results["log2FoldChange"].abs() >= lfc)
    up = significant & (results["log2FoldChange"] > 0)
    down = significant & (results["log2FoldChange"] < 0)
    write_bed(results, up, f"{output_prefix}_{treatment_condition}.bed")
    write_bed(results, down, f"{output_prefix}_{reference_condition}.bed")
    plot_results(results, count_matrix, factors, samples, output_prefix, fdr, lfc)

    diagnostics = pd.DataFrame(count_diagnostics).merge(bw_diagnostics, on="sample")
    diagnostics.to_csv(diagnostic_path, sep="\t", index=False)
    metadata = {
        "method": "DESeq2-style negative-binomial Wald test with PAW-derived region/sample offsets",
        "regions": str(Path(regions_file).resolve()),
        "sample_sheet": str(Path(sample_sheet).resolve()),
        "reference_condition": reference_condition,
        "treatment_condition": treatment_condition,
        "mapq": mapq,
        "min_count": min_count,
        "fdr": fdr,
        "absolute_log2_fold_change": lfc,
        "tested_regions": int(np.isfinite(results["pvalue"]).sum()),
        "significant_reference_enriched": int(down.sum()),
        "significant_treatment_enriched": int(up.sum()),
        "runtime": str(datetime.now() - start_time),
        "warning": "This is a DESeq2-style implementation, not exact DESeq2/PyDESeq2 output.",
    }
    with open(f"{output_prefix}_run.json", "w") as handle:
        json.dump(metadata, handle, indent=2)
    rprint(f"[{output_prefix}] finished: {down.sum():,} {reference_condition}-enriched; {up.sum():,} {treatment_condition}-enriched")


if __name__ == "__main__":
    patrol_replicates()
