#!/usr/bin/env python
# --coding:utf-8--

"""
PATROL algorithm for cross-sample normalized ChIP-seq/DNase-seq/ATAC-seq differential peaks detection.

This tool detects highly variable features by performing a two-pass Mahalanobis distance test on normalized 
signal data. The algorithm uses internal references from DNase-seq, ATAC-seq, ChIP-seq, and MNase-seq experiments.

Author: CAO Yaqiang
Date: 2025-04-09
Email: caoyaqiang0410@gmail.com
"""

#sys
import os
import warnings
from datetime import datetime
warnings.filterwarnings("ignore")

#3rd
import click
import pyBigWig
import numpy as np
import pandas as pd
from tqdm import tqdm
from scipy.stats import chi2
from joblib import Parallel, delayed

# Plotting libraries and settings
import matplotlib as mpl
mpl.use("pdf")
import seaborn as sns
import pylab
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["figure.figsize"] = (3.2, 2.2)
mpl.rcParams["savefig.transparent"] = True
mpl.rcParams["savefig.bbox"] = "tight"
mpl.rcParams["font.size"] = 7.0
mpl.rcParams["font.sans-serif"] = "Arial"
mpl.rcParams["savefig.format"] = "pdf"
sns.set_style("white")
# Define a list of colors for plotting
colors = [
    (0.8941176470588236, 0.10196078431372549, 0.10980392156862745),
    (0.21568627450980393, 0.49411764705882355, 0.7215686274509804),
    (0.30196078431372547, 0.6862745098039216, 0.2901960784313726),
    (0.596078431372549, 0.3058823529411765, 0.6392156862745098),
    (1.0, 0.4980392156862745, 0.0),
    (0, 0, 0),  # black
    (0.6509803921568628, 0.33725490196078434, 0.1568627450980392),
    (0.9686274509803922, 0.5058823529411764, 0.7490196078431373),
    (0.6, 0.6, 0.6),
    (0.5529411764705883, 0.8274509803921568, 0.7803921568627451),
    (1.0, 1.0, 0.7019607843137254),
    (0.7450980392156863, 0.7294117647058823, 0.8549019607843137),
    (0.984313725490196, 0.5019607843137255, 0.4470588235294118),
    (0.5019607843137255, 0.6941176470588235, 0.8274509803921568),
    (0.9921568627450981, 0.7058823529411765, 0.3843137254901961),
    (0.7019607843137254, 0.8705882352941177, 0.4117647058823529),
    (0.9882352941176471, 0.803921568627451, 0.8980392156862745),
]


######################################
# Utility Functions                  #
######################################
def rprint(message):
    """
    Print a message with the current datetime prepended.

    :param message: str, the message to print.
    :return: None
    """
    report = "\t".join([str(datetime.now()), message])
    print(report)


def readBed(filepath):
    """
    Read a BED format file and return regions with at least three columns.
    Regions with a width less than 100 bp are ignored.

    :param filepath: str, path to the BED file.
    :return: list of regions, each as [chrom, start, end].
    """
    regions = []
    with open(filepath) as f:
        for line in f:
            fields = line.strip().split("\t")
            if len(fields) < 3:
                continue
            try:
                fields[1] = int(fields[1])
                fields[2] = int(fields[2])
            except Exception:
                continue
            if fields[2] - fields[1] < 100:
                continue
            regions.append(fields[:3])
    return regions


######################################
# Region Handling Functions          #
######################################
def buildCov(regions):
    """
    Build a coverage dictionary and range limits for genomic regions.

    :param regions: list of regions [chrom, start, end].
    :return: tuple (cov, lims)
             - cov: dict mapping each chromosome to a set of covered positions.
             - lims: dict mapping each chromosome to [min_start, max_end].
    """
    cov = {}
    lims = {}
    for region in regions:
        chrom, start, end = region
        cov.setdefault(chrom, set()).update(range(start, end + 1))
        if chrom not in lims:
            lims[chrom] = [start, end]
        else:
            if start < lims[chrom][0]:
                lims[chrom][0] = start
            if end > lims[chrom][1]:
                lims[chrom][1] = end
    return cov, lims


def getRegion(cov, margin=1):
    """
    Merge individual covered positions into continuous regions.

    :param cov: dict mapping chromosome to a set of covered positions.
    :param margin: int, maximum gap allowed between positions to be merged.
    :return: list of regions as [chrom, start, end].
    """
    regions = []
    for chrom, pos_set in cov.items():
        pos_list = sorted(list(pos_set))
        i = 0
        while i < len(pos_list) - 1:
            # Find the end of a continuous region
            for j in range(i + 1, len(pos_list)):
                if pos_list[j] - pos_list[j - 1] > margin:
                    break
            start = pos_list[i]
            end = pos_list[j - 1]
            regions.append([chrom, start, end])
            i = j
    return regions


def checkBgOverlaps(chrom, start, end, cov, lims):
    """
    Check if a proposed background region overlaps with reference coverage or falls outside limits.

    :param chrom: str, chromosome name.
    :param start: int, start coordinate.
    :param end: int, end coordinate.
    :param cov: dict mapping chromosome to a set of covered positions.
    :param lims: dict mapping chromosome to [min_start, max_end].
    :return: bool, True if overlap exists or region is out-of-bound.
    """
    if start < lims[chrom][0] or end > lims[chrom][1]:
        return True
    for pos in range(start, end + 1):
        if pos in cov[chrom]:
            return True
    return False


def getFgBgs(refPeaks, exts=[10]):
    """
    Generate foreground (signal) regions and corresponding background regions.

    :param refPeaks: list of reference peaks [chrom, start, end].
    :param exts: list of int, extension multipliers to generate background regions.
    :return: tuple (foreground_regions, background_regions).
             Foreground regions are the merged reference regions,
             while background regions are generated by shifting the foreground regions.
    """
    refCov, refLims = buildCov(refPeaks)
    fg_regions = getRegion(refCov)
    bg_regions = []
    for region in fg_regions:
        chrom, start, end = region
        d = end - start
        for ext in exts:
            # Extend in the negative direction
            s = start - ext * d
            e = end - ext * d
            if checkBgOverlaps(chrom, s, e, refCov, refLims):
                continue
            bg_regions.append([chrom, s, e])
            # Extend in the positive direction
            s = start + ext * d
            e = end + ext * d
            if checkBgOverlaps(chrom, s, e, refCov, refLims):
                continue
            bg_regions.append([chrom, s, e])
    return fg_regions, bg_regions


######################################
# Quantification Functions           #
######################################
def quant(regions, bg_regions, bw_filepath):
    """
    Quantify genomic signals from a bigWig file over specified regions.
    Background regions are used to estimate a noise level to avoid division by zero.

    :param regions: list of regions [chrom, start, end] of interest.
    :param bg_regions: list of background regions [chrom, start, end].
    :param bw_filepath: str, path to the bigWig file.
    :return: tuple (signal_series, noise_level)
             - signal_series: pandas Series mapping region identifier to averaged signal.
             - noise_level: float, average noise signal computed from background regions.
    """
    bw = pyBigWig.open(bw_filepath)
    noise_vals = []
    # Calculate noise level from background regions
    for region in bg_regions:
        chrom, start, end = region
        vals = bw.values(chrom, start, end)
        vals = np.nan_to_num(vals)
        noise_vals.append(np.sum(vals) / len(vals))
    noise_level = np.mean(noise_vals)

    signals = {}
    for region in tqdm(regions):
        chrom, start, end = region
        region_id = f"{chrom}:{start}-{end}"
        vals = bw.values(chrom, start, end)
        vals = np.nan_to_num(vals)
        avg_signal = np.sum(vals) / len(vals)
        # Replace zero values with the estimated noise level
        if avg_signal == 0:
            avg_signal = noise_level
        signals[region_id] = avg_signal
    bw.close()
    return pd.Series(signals), noise_level


######################################
# Statistical Functions              #
######################################
def mahalanobis(matrix):
    """
    Calculate the Mahalanobis distance for each row of a matrix and compute corresponding chi-square p-values.

    :param matrix: numpy array where rows are observations and columns are variables.
    :return: tuple (distances, p_values)
             - distances: numpy array of Mahalanobis distances.
             - p_values: numpy array of chi-square p-values.
    """
    cov = np.cov(matrix, rowvar=False)
    invCov = np.linalg.inv(cov)
    center = np.mean(matrix, axis=0)
    diff = matrix - center
    distances = np.dot(np.dot(diff, invCov), diff.T).diagonal()
    p_values = 1 - chi2.cdf(distances, matrix.shape[1] - 1)
    return distances, p_values


def twoPassesMDTest(data, pcut=0.01):
    """
    Perform a two-pass Mahalanobis distance test to detect outliers.

    First pass uses all data to flag outliers.
    Second pass recomputes the covariance matrix and center using only non-outlier data.

    :param data: pandas DataFrame where rows are observations and columns are features.
    :param pcut: float, p-value cutoff for flagging outliers.
    :return: tuple (distances, p_values) for all data based on recomputed statistics.
    """
    # First pass: calculate distances and p-values
    distances, p_values = mahalanobis(data.values)
    distances = pd.Series(distances, index=data.index)
    p_values = pd.Series(p_values, index=data.index)
    # Identify outliers
    outlier_indices = p_values[p_values < pcut].index
    # Second pass: recalculate using data without outliers
    refined_data = data.drop(outlier_indices)
    cov = np.cov(refined_data.values, rowvar=False)
    invCov = np.linalg.inv(cov)
    center = np.mean(refined_data.values, axis=0)
    diff = data.values - center
    distances = np.dot(np.dot(diff, invCov), diff.T).diagonal()
    p_values = 1 - chi2.cdf(distances, data.shape[1] - 1)
    distances = pd.Series(distances, index=data.index)
    p_values = pd.Series(p_values, index=data.index)
    return distances, p_values


######################################
# Plotting Functions                 #
######################################
def showMA(m, a, control_label, treatment_label, p_values, out_prefix, noise, pcut=0.01):
    """
    Generate and save an MA plot showing differential signal (M) versus average signal (A).

    Points are colored based on whether they pass the specified p-value cutoff.

    :param m: pandas Series, log2 fold-change (treatment - control).
    :param a: pandas Series, average log2 signal.
    :param control_label: str, label for control sample.
    :param treatment_label: str, label for treatment sample.
    :param p_values: pandas Series, p-values from the Mahalanobis test.
    :param out_prefix: str, prefix for the output plot file.
    :param noise: float, estimated noise level for scaling.
    :param pcut: float, p-value cutoff to highlight significant regions.
    """
    fig, ax = pylab.subplots()
    # Plot all regions in gray
    ax.scatter(a, m, s=0.5, color="gray", alpha=0.6, label=f"Total {len(a)} regions")
    
    # Identify regions passing the cutoff and above noise
    significant = p_values[p_values <= pcut].index
    valid = a[significant][a[significant] > noise].index
    m_valid = m[valid]
    # Split based on sign of M (up- or down-regulated)
    up_idx = m_valid[m_valid > 0].index
    down_idx = m_valid[m_valid < 0].index
    ax.scatter(a[down_idx], m[down_idx], s=1, color=colors[0], alpha=0.8, label=f"{len(down_idx)} {control_label} regions")
    ax.scatter(a[up_idx], m[up_idx], s=1, color=colors[1], alpha=0.8, label=f"{len(up_idx)} {treatment_label} regions")
    leg = ax.legend(loc='best')
    for handle in leg.legendHandles:
        handle._sizes = [10]
    ax.axhline(0, color="gray", linestyle="--")
    ax.set_title(f"Signal Comparison\nMahalanobis P-value < {pcut}")
    ax.set_xlabel(f"log2({treatment_label}) + log2({control_label})")
    ax.set_ylabel(f"log2({treatment_label}) - log2({control_label})")
    pylab.savefig(out_prefix + "_MA.pdf")


######################################
# Main Workflow (CLI)                #
######################################
@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.option(
    "-r",
    required=True,
    type=str,
    help=(
        "Path to a BED format file specifying regions of interest for differential analysis. "
        "Each region must have at least three columns: chromosome, start, and end. "
        "Example: regions.bed"
    ),
)
@click.option(
    "-o",
    required=True,
    type=str,
    help=(
        "Output file prefix. Several output files (statistics and plots) will be generated using this prefix. "
        "Example: results/test"
    ),
)
@click.option(
    "-c",
    required=True,
    type=str,
    help=(
        "Path to the control sample bigWig file. This file should be normalized (e.g., to RPM). "
        "Example: control.bw"
    ),
)
@click.option(
    "-lc",
    default="control",
    type=str,
    help="Label for the control sample (default: 'control').",
)
@click.option(
    "-t",
    required=True,
    type=str,
    help=(
        "Path to the treatment sample bigWig file. This file should be normalized (e.g., to RPM). "
        "Example: treatment.bw"
    ),
)
@click.option(
    "-lt",
    default="trt",
    type=str,
    help="Label for the treatment sample (default: 'trt').",
)
@click.option(
    "-pcut",
    default=0.05,
    type=float,
    help=(
        "P-value cutoff for the two-pass Mahalanobis distance test. Regions with p-values below this cutoff "
        "are considered outliers. Default is 0.05."
    ),
)
def patrol(r, c, t, o, lc, lt, pcut):
    """
    PATROL: Sequencing Data Variable Features Detection Algorithm.

    This tool detects highly variable genomic features by performing a two-pass Mahalanobis distance test on normalized ChIP-seq/DNase-seq/ATAC-seq data. It quantifies signals over specified regions, estimates background noise, and identifies differential signals using statistical tests.

    Examples:
      $ patrol.py -r regions.bed -c control.bw -t treatment.bw -o results_diff

    """
    start_time = datetime.now()
    script_name = os.path.basename(__file__)
    rprint(f"{script_name} -r {r} -o {o} -c {c} -t {t} -lc {lc} -lt {lt} -pcut {pcut}")

    # Step 0: Check required input files exist
    for filepath in [r, c, t]:
        if not os.path.isfile(filepath):
            rprint(f"ERROR! Input file {filepath} does not exist. Exiting.")
            return
    output_dir = os.path.dirname(o)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)

    # Step 1: Read regions of interest and generate background regions to estimate noise
    regions = readBed(r)
    fg_regions, bg_regions = getFgBgs(regions)
    
    # Step 2: Quantify signals for control and treatment samples
    rprint(f"[{o}] Step 1: Quantifying signal for control sample ({lc})")
    control_signals, control_noise = quant(fg_regions, bg_regions, c)
    rprint(f"[{o}] Step 1: Quantifying signal for treatment sample ({lt})")
    treatment_signals, treatment_noise = quant(fg_regions, bg_regions, t)
    
    # Step 3: Compute log2-transformed signals and derive M (difference) and A (average) values
    control_log = np.log2(control_signals)
    treatment_log = np.log2(treatment_signals)
    A = (control_log + treatment_log) / 2
    M = treatment_log - control_log
    md_data = pd.DataFrame({"m": M, "a": A})
    
    # Step 4: Perform two-pass Mahalanobis distance test to detect outliers
    rprint(f"[{o}] Step 2: Performing two-pass Mahalanobis distance test with pcut = {pcut}")
    distances, p_values = twoPassesMDTest(md_data, pcut=pcut)
    
    # Save statistics to file
    stats_data = {
        lc: control_signals, 
        lt: treatment_signals, 
        f"log2 fold change ({lt}/{lc})": M, 
        "average": A,
        "MD": distances, 
        "Chi-squared P-value": p_values
    }
    stats_df = pd.DataFrame(stats_data)
    stats_df.to_csv(f"{o}_stat.txt", sep="\t", index_label="regionId")
    rprint(f"[{o}] Step 2: Statistical results saved to {o}_stat.txt.")

    # Step 5: Generate and save MA plot of differential signals
    rprint(f"[{o}] Step 3: Generating MA plot of detected differential features")
    combined_noise = A[p_values[p_values>pcut].index].min()
    showMA(M, A, lc, lt, p_values, o, combined_noise, pcut=pcut)

    # Step 6: Show aggregated differential peaks 

    end_time = datetime.now()
    rprint(f"{script_name} job finished. Used time: {end_time - start_time}")


if __name__ == "__main__":
    patrol()
