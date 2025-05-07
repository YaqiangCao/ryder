#!/usr/bin/env python
# --coding:utf-8--
"""
PAW algorithm for cross-sample ChIP-seq/DNase-seq/ATAC-seq normalization with internal control. 

Normalizes ChIP-seq, DNase-seq, or ATAC-seq data between a target sample and a control/reference sample using internal reference regions (e.g., housekeeping gene TSSs, conserved CTCF sites) assumed to be stable across conditions.

The algorithm involves:
1. Identifying stable reference regions by removing outliers using Mahalanobis distance.
2. Defining foreground (reference) and background regions.
3. Calculating initial QC metrics (signal-to-noise, background levels).
4. Estimating scaling factors for background noise and reference signal region.
5. Estimating linear transformation parameters (alpha, beta) for signal regions after log2 transformation.
6. Training a Gaussian Mixture Model (GMM) to distinguish signal from noise in the target sample based on reference/background regions.
7. Applying normalization to the target BigWig file:
    - Noise regions are scaled by the background factor.
    - Signal regions are transformed using the estimated alpha and beta on log2-scaled data with reference scaling factor, then converted back.

Assumptions:
    - Input BigWig files are pre-normalized (e.g., to RPM, CPM).
    - Provided reference regions contain a sufficient number of stable sites.
    - Background and signal log2 distributions are normal and have same std.
History:
2025-04-04: Initial version
2025-04-09: Revised and basically finished
2025-04-11: Introduced Mahalanobis distance test for outlier removal in reference regions. 
2025-04-30: Add -pred option for bgs, alpha, beta pretrained from spike-in data
2025-05-01: Add -mode option for linear fitting or normal distribution draw
2025-05-02: Integrate signal region scaling factor
2025-05-02: remove GMM and estimate a noise cutoff, much faster and general applied
"""

__author__ = "CAO Yaqiang"
__date__ = "2025-04-04"
__modified__ = ""
__email__ = "caoyaqiang0410@gmail.com"

#sys
import os
import random
import warnings
from glob import glob
from datetime import datetime

warnings.filterwarnings("ignore")

#3rd
import click
import pyBigWig
import numpy as np
import pandas as pd
from tqdm import tqdm
from scipy.stats import chi2
from sklearn import linear_model
from joblib import Parallel, delayed
from sklearn.mixture import GaussianMixture as GMM
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_absolute_error as MAE

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
    Read a BED format file and return regions with at least 3 columns.
    Regions with width less than 100 are skipped.
    
    :param filepath: str, path to BED file.
    :return: list of regions, each region is a list: [chrom, start, end]
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
# I/O and Quantification Functions   #
######################################
def quant(regions, bw_filepath):
    """
    Quantify genomic features from a bigWig file for given regions.

    :param regions: list of regions; each region is a list [chrom, start, end].
    :param bw_filepath: str, path to bigWig file.
    :return: pandas Series mapping region identifier ("chrom:start-end") to averaged signal value.
    """
    bw = pyBigWig.open(bw_filepath)
    result = {}
    for region in tqdm(regions):
        chrom, start, end = region
        region_id = f"{chrom}:{start}-{end}"
        values = bw.values(chrom, start, end)
        # Replace NaNs with zeros and compute mean signal
        values = np.nan_to_num(values)
        result[region_id] = np.sum(values) / len(values)
    bw.close()
    return pd.Series(result)


def getBinMean(signal_array, bins=100):
    """
    Compute mean signal values over bins for a given array.

    :param signal_array: numpy array of signal values.
    :param bins: int, number of bins.
    :return: numpy array of bin means.
    """
    width = int(len(signal_array) / bins)
    binned_mean = signal_array[:bins * width].reshape(-1, width).mean(axis=1)
    return binned_mean


def _getBwSig(region, bw_filepath, bins=100):
    """
    Fetch signal values from a bigWig file for a given region and return the binned mean signal.
    
    :param region: list, [chrom, start, end].
    :param bw_filepath: str, path to bigWig file.
    :param bins: int, number of bins.
    :return: numpy array of binned signal means or zeros if an error occurs.
    """
    bw = pyBigWig.open(bw_filepath)
    try:
        values = bw.values(region[0], region[1], region[2])
        values = np.nan_to_num(values)
        if len(values) > bins:
            values = getBinMean(values, bins)
            return values
        else:
            return np.zeros(bins)
    except Exception:
        return np.zeros(bins)
    finally:
        bw.close()


def getBwSig(regions, bw_filepath, bins=100, n_jobs=2):
    """
    Get binned signal values for multiple regions from a bigWig file in parallel.

    :param regions: list of regions; each region is [chrom, start, end].
    :param bw_filepath: str, path to bigWig file.
    :param bins: int, number of bins for conversion.
    :param n_jobs: int, number of CPUs for parallel processing.
    :return: numpy array of signals.
    """
    signals = Parallel(n_jobs=n_jobs, backend="multiprocessing")(
        delayed(_getBwSig)(region, bw_filepath, bins) for region in regions)
    return np.array(signals)


######################################
# Statistical Functions              #
######################################
def mahalanobis(matrix):
    """
    Calculate the Mahalanobis distance for each row in the matrix and corresponding chi-square test p-values.
    
    :param matrix: numpy array, each row represents an item and each column a variable.
    :return: tuple (dis, ps)
             - dis: numpy array, Mahalanobis distances.
             - ps: numpy array, chi-square test p-values.
    """
    cov = np.cov(matrix, rowvar=False)
    invCov = np.linalg.inv(cov)
    center = np.mean(matrix, axis=0)
    diff = matrix - center
    distances = np.dot(np.dot(diff, invCov), diff.T).diagonal()
    p_values = 1 - chi2.cdf(distances, matrix.shape[1] - 1)
    return distances, p_values


######################################
# Plotting and Visualization         #
######################################
def showSig(regions,
            control_bw,
            treatment_bw,
            out_prefix,
            title="",
            refLabel="ref",
            tgtLabel="tgt",
            ext=10000):
    """
    Plot and save the averaged signal around the center of regions for both control and treatment.

    :param regions: list of regions [chrom, start, end].
    :param control_bw: str, bigWig file path for the control sample.
    :param treatment_bw: str, bigWig file path for the treatment sample.
    :param out_prefix: str, output file prefix.
    :param title: str, plot title.
    :param refLabel: str, label for control sample.
    :param tgtLabel: str, label for treatment sample.
    :param ext: int, extension (in bp) around the center.
    """
    x_vals = np.arange(-ext, ext + 1)
    control_signal = np.zeros(len(x_vals))
    treatment_signal = np.zeros(len(x_vals))
    bw_control = pyBigWig.open(control_bw)
    bw_treatment = pyBigWig.open(treatment_bw)
    for region in tqdm(regions):
        chrom = region[0]
        center = int((region[1] + region[2]) / 2)
        control_vals = bw_control.values(chrom, center - ext, center + ext + 1)
        control_vals = np.nan_to_num(control_vals)
        control_signal += control_vals
        treatment_vals = bw_treatment.values(chrom, center - ext,
                                             center + ext + 1)
        treatment_vals = np.nan_to_num(treatment_vals)
        treatment_signal += treatment_vals
    control_signal /= len(regions)
    treatment_signal /= len(regions)
    fig, ax = pylab.subplots()
    ax.plot(x_vals, control_signal, label=refLabel, color=colors[0])
    ax.plot(x_vals, treatment_signal, label=tgtLabel, color=colors[1])
    ax.set_ylabel("avg. signal")
    ax.set_xlabel("distance from reference center")
    ax.set_title(title)
    ax.set_xlim([-int(ext / 2), int(ext / 2) + 1])
    ax.legend(frameon=False, markerscale=0, labelcolor=[colors[0], colors[1]])
    pylab.savefig(out_prefix + ".pdf")


def plotBinSig(ax, x, y1, y2, label1, label2, xlabel, ylabel, title):
    """
    Plot two binned signal lines on a given axis.

    :param ax: matplotlib axis object.
    :param x: x-axis values.
    :param y1: y-axis values for the first signal.
    :param y2: y-axis values for the second signal.
    :param label1: label for first signal.
    :param label2: label for second signal.
    :param xlabel: label for the x-axis.
    :param ylabel: label for the y-axis.
    :param title: plot title.
    :return: the axis object with the plot.
    """
    ax.plot(x, y1, label=label1, color=colors[0])
    ax.plot(x, y2, label=label2, color=colors[1])
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend(frameon=False, markerscale=0, labelcolor=[colors[0], colors[1]])
    return ax


def getQc(fg_regions,
          bg_regions,
          ref_bw,
          tgt_bw,
          out_prefix,
          refLabel="ref",
          tgtLabel="tgt",
          bins=100):
    """
    Quality control: compare noise and signal-to-noise ratio between reference and target samples.

    :param fg_regions: list of foreground (reference) regions.
    :param bg_regions: list of background regions.
    :param ref_bw: str, bigWig file for the reference sample.
    :param tgt_bw: str, bigWig file for the target sample.
    :param out_prefix: str, output file prefix.
    :param refLabel: str, label for reference.
    :param tgtLabel: str, label for target.
    :param bins: int, number of bins.
    :return: tuple (fgRef, fgTgt, bgRef, bgTgt) each a numpy array of averaged binned signals.
    """
    fgRef = getBwSig(fg_regions, ref_bw, bins=bins)
    fgRef = fgRef.mean(axis=0)
    fgTgt = getBwSig(fg_regions, tgt_bw, bins=bins)
    fgTgt = fgTgt.mean(axis=0)
    bgRef = getBwSig(bg_regions, ref_bw, bins=bins)
    bgRef = bgRef.mean(axis=0)
    bgTgt = getBwSig(bg_regions, tgt_bw, bins=bins)
    bgTgt = bgTgt.mean(axis=0)

    bg_scale = bgRef.mean() / bgTgt.mean()
    fg_scale = fgRef.mean() / fgTgt.mean()

    fig, axs = pylab.subplots(1, 3, figsize=(6, 2), sharex=True)
    axs = axs.reshape(-1)
    x = np.arange(bins)

    # Background level plot
    ax = axs[0]
    plotBinSig(ax, x, bgRef, bgTgt, refLabel, tgtLabel, "bins", "avg. signal",
               f"background\n sf:%.3f" % bg_scale)

    # Signal level plot
    ax = axs[1]
    plotBinSig(ax, x, fgRef, fgTgt, refLabel, tgtLabel, "bins", "avg. signal",
               f"reference region\n sf:%.3f" % fg_scale)

    # Signal-to-noise ratio plot
    refSN = fgRef / bgRef
    tgtSN = fgTgt / bgTgt
    ax = axs[2]
    sn_scale = refSN.mean() / tgtSN.mean()
    plotBinSig(ax, x, refSN, tgtSN, refLabel, tgtLabel, "bins", "S2N ratio",
               f"signal to noise ratio\n sf:%.3f" % sn_scale)

    pylab.tight_layout()
    pylab.savefig(out_prefix + "_2_qc.pdf")
    return fgRef, fgTgt, bgRef, bgTgt


def estLAB(refs, tgts, cpu=1):
    """
    Estimate the linear regression alpha and beta of target to control.

    :param refs: numpy.array, reference signal
    :param tgts: numpy.array, target signal 
    :param cpu: int, cpu number to run the fitting jobs
    """
    x = [[t] for t in tgts]
    x = np.array(x)
    y = np.array(refs)
    x_train, x_vali, y_train, y_vali = train_test_split(x, y, test_size=0.1)
    lra = linear_model.LinearRegression(n_jobs=cpu, fit_intercept=True)
    lrb = linear_model.TheilSenRegressor(random_state=123,
                                         n_jobs=cpu,
                                         fit_intercept=True)
    lra.fit(x_train, y_train)
    lrb.fit(x_train, y_train)
    xv = np.array([t[0] for t in x_vali])
    ypa = lra.predict(x_vali)
    ypb = lrb.predict(x_vali)
    maea = MAE(y_vali, ypa)
    maeb = MAE(y_vali, ypb)
    if maea < maeb:
        lr = lra
    else:
        lr = lrb
    alpha = lr.coef_[0]
    beta = lr.intercept_
    return alpha, beta


def estFit(fg_regions,
           ref_bw,
           tgt_bw,
           out_prefix,
           fg_sf,
           refLabel="ref",
           tgtLabel="tgt",
           alpha=None,
           beta=None,
           mode="norm",
           cpu=1,
           bins=5):
    """
    Estimate the linear fit parameters between the reference and target sample signal regions.

    :param fg_regions: list of foreground regions.
    :param ref_bw: str, bigWig file for the reference sample.
    :param tgt_bw: str, bigWig file for the target sample.
    :param out_prefix: str, output file prefix.
    :param refLabel: str, label for reference.
    :param tgtLabel: str, label for target.
    :param fg_sf: float, scaling factor with raw data
    :param alpha: float, linear fitting alpha, if None, estimate it.
    :param beta: float, linear fitting beta, if None, estimate it.
    :param mode: str, norm or lr, norm means normal distribution draw and lr means linear fitting
    :param bins: int, number of bins.
    :return: list [alpha, beta] for linear fit: log2(target) = alpha * log2(target) + beta.
    """
    refS = getBwSig(fg_regions, ref_bw, bins=bins)
    refS = pd.Series(refS.reshape(-1))
    tgtS = getBwSig(fg_regions, tgt_bw, bins=bins)
    tgtS = pd.Series(tgtS.reshape(-1))

    fig, axs = pylab.subplots(2, 2, figsize=(4, 4))
    axs = axs.reshape(-1)

    # Filter out zeros and perform log2 transformation
    idx = refS[refS > 0].index
    idx = tgtS[idx][tgtS[idx] > 0].index
    refS = refS[idx]
    tgtS = tgtS[idx]
    refS = np.log2(refS)
    tgtS = np.log2(tgtS)

    # Plot signal distributions
    ax = axs[0]
    sns.kdeplot(refS, label=refLabel, ax=ax, fill=True, color=colors[0])
    sns.kdeplot(tgtS, label=tgtLabel, ax=ax, fill=True, color=colors[1])
    ax.legend(frameon=False, markerscale=0, labelcolor=[colors[0], colors[1]])
    ax.set_xlabel("log2(signal)")
    ax.set_title("signal distribution")

    # Plot linear fitting before correction
    ax = axs[1]
    ax.scatter(refS, tgtS, s=0.1)
    diff = tgtS - refS
    avg = (refS + tgtS) / 2
    pcc = diff.corr(avg)
    s_min = np.min([refS.min(), tgtS.min()])
    s_max = np.max([refS.max(), tgtS.max()])
    ax.plot([s_min, s_max], [s_min, s_max], color="gray", linestyle="--")
    ax.set_title(f"before normalization \n M~A PCC:{pcc:.3f}")
    ax.set_xlabel(refLabel)
    ax.set_ylabel(tgtLabel)

    #just apply a global scaling factor
    tgt_scaled = np.log2((2**tgtS) * fg_sf)
    if alpha is None and beta is None:
        if mode == "norm":
            #assuming normal distribution and just draw the distribution to same
            alpha = refS.std() / tgt_scaled.std()
            beta = refS.mean() - alpha * tgt_scaled.mean()
        if mode == "lr":
            alpha, beta = estLAB(refS, tgt_scaled, cpu=cpu)
    tgt_scaled = tgt_scaled * alpha + beta
    diff_scaled = tgt_scaled - refS
    avg_scaled = (refS + tgt_scaled) / 2

    # Plot signal distribution after correction
    ax = axs[2]
    sns.kdeplot(refS, label=refLabel, ax=ax, fill=True, color=colors[0])
    sns.kdeplot(tgt_scaled, label=tgtLabel, ax=ax, fill=True, color=colors[1])
    ax.legend(frameon=False, markerscale=0, labelcolor=[colors[0], colors[1]])
    ax.set_xlabel("log2(signal)")
    ax.set_title("after normalization")

    # Distribution matching (signal conversion)
    ax = axs[3]
    ax.scatter(avg_scaled, diff_scaled, s=0.1)
    if beta > 0:
        ax.set_title(
            "after normalization M~A PCC:%.3f\nlog2(%s)=%.3flog2(%.3f%s)+%.3f"
            % (diff_scaled.corr(avg_scaled), tgtLabel, alpha, fg_sf, tgtLabel,
               beta))
    else:
        ax.set_title(
            "after normalization M~A PCC:%.3f\nlog2(%s)=%.3flog2(%.3f%s)%.3f" %
            (diff_scaled.corr(avg_scaled), tgtLabel, alpha, fg_sf, tgtLabel,
             beta))
    ax.axhline(0, color="gray", linestyle="--")
    ax.set_xlabel(f"A, (log2({refLabel})+log2({tgtLabel}))/2")
    ax.set_ylabel(f"M, log2({tgtLabel})-log2({refLabel})")

    pylab.tight_layout()
    pylab.savefig(out_prefix + "_3_fgSignalConversion.pdf")
    return [alpha, beta]



def getNoiseCut(fg_regions, bg_regions, bw_filepath, out_prefix, ext=5000, bins=5):
    """
    Get a noise level cutoff to classify background and signal regions for the target sample.

    :param fg_regions: list of foreground regions.
    :param bg_regions: list of background regions.
    :param tgt_bw: str, path to target sample bigWig file.
    :param out_prefix: str, output file prefix.
    :param tgtLabel: str, label for target sample.
    :param ext: int, extension size for nearby regions (bp).
    :return: tuple (tgtGmm, tgt_class_mapping)
    """
 
    # Get signals and apply log2 transformation on positive values only
    fg_signal = getBwSig(fg_regions, bw_filepath, bins=bins).reshape(-1)
    fg_signal = np.log2(fg_signal[fg_signal > 0])
    bg_signal = getBwSig(bg_regions, bw_filepath, bins=bins).reshape(-1)
    bg_signal = np.log2(bg_signal[bg_signal > 0])
    #assuming two normal distribution and std is same and the intersection is the mean of two
    noise = (np.mean(fg_signal) + np.mean(bg_signal)) / 2
    fig, ax = pylab.subplots()
    sns.kdeplot(bg_signal, label="background", fill=False, ax=ax)
    sns.kdeplot(fg_signal, label="reference region", fill=False, ax=ax)
    ax.axvline(x=noise, color="gray", linestyle="--")
    ax.legend(frameon=False, markerscale=0)
    ax.set_xlabel("log2(signal)")
    pylab.savefig(out_prefix + "_4_NoiseCutoff.pdf")
    return 2**noise


######################################
# Region Handling Functions          #
######################################
def buildCov(regions):
    """
    Build a coverage dictionary and range limitations for genomic regions.
    
    :param regions: list of regions [chrom, start, end].
    :return: tuple (cov, lims)
             - cov: dict mapping chromosome to a set of covered positions.
             - lims: dict mapping chromosome to [min_start, max_end] limits.
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


def getRegion(coverage):
    """
    Get continuous regions from a coverage dictionary.

    :param coverage: dict mapping chromosome to a set of covered positions.
    :return: list of regions [chrom, start, end].
    """
    regions = []
    for chrom, pos_set in coverage.items():
        pos_list = sorted(list(pos_set))
        i = 0
        while i < len(pos_list) - 1:
            for j in range(i + 1, len(pos_list)):
                if pos_list[j] - pos_list[j - 1] > 1:
                    break
            start = pos_list[i]
            end = pos_list[j - 1]
            i = j
            regions.append([chrom, start, end])
    return regions


def checkBgOverlaps(chrom, start, end, coverage, limits):
    """
    Check if a background region overlaps with the reference coverage or exceeds defined limits.

    :param chrom: str, chromosome name.
    :param start: int, start coordinate.
    :param end: int, end coordinate.
    :param coverage: dict mapping chromosome to a set of covered positions.
    :param limits: dict mapping chromosome to [min_start, max_end].
    :return: bool, True if there is an overlap or region is out-of-bound.
    """
    if start < limits[chrom][0] or end > limits[chrom][1]:
        return True
    for pos in range(start, end + 1):
        if pos in coverage[chrom]:
            return True
    return False


def getFgBgs(refPeaks):
    """
    Generate foreground (signal) and background regions based on reference peaks.
    
    :param refPeaks: list of reference peaks (regions) [chrom, start, end].
    :return: tuple (fg_regions, bg_regions)
             - fg_regions: list of sorted reference regions.
             - bg_regions: list of background regions generated by extending fg_regions.
    """
    refCov, refLims = buildCov(refPeaks)
    fg_regions = getRegion(refCov)
    bg_regions = []
    for region in fg_regions:
        chrom, start, end = region
        d = end - start
        for ext in [1, 2, 3,4, 5, 10]:
            s = start - ext * d
            e = end - ext * d
            if checkBgOverlaps(chrom, s, e, refCov, refLims):
                continue
            bg_regions.append([chrom, s, e])
            s = start + ext * d
            e = end + ext * d
            if checkBgOverlaps(chrom, s, e, refCov, refLims):
                continue
            bg_regions.append([chrom, s, e])
    return fg_regions, bg_regions


######################################
# Outlier Removal Function           #
######################################
def removeOutliers(refPeaks, control_bw, treatment_bw, pcut=0.1):
    """
    Filter reference peaks by removing outliers using the Mahalanobis distance test.
    
    Steps:
      1. Quantify signal for control and treatment samples over the reference peaks.
      2. Filter out regions with non-positive signal.
      3. Compute log2-transformed signals and derive average (A) and difference (M).
      4. Calculate Mahalanobis distances and corresponding p-values.
      5. Retain regions with p-value > pcut.
    
    :param refPeaks: list of reference peaks.
    :param control_bw: str, path to control bigWig file.
    :param treatment_bw: str, path to treatment bigWig file.
    :param pcut: float, p-value cutoff for filtering.
    :return: list of filtered regions.
    """
    # Quantify signals from control and treatment bigWig files
    sc = quant(refPeaks, control_bw)
    st = quant(refPeaks, treatment_bw)
    sc = sc[sc > 0]
    st = st[st > 0]
    common_index = sc.index.intersection(st.index)
    sc = np.log2(sc[common_index])
    st = np.log2(st[common_index])
    # Construct dataframe with average (A) and difference (M)
    data = pd.DataFrame({"a": (sc + st) / 2, "m": (st - sc)})
    # Compute Mahalanobis distances and p-values
    distances, p_values = mahalanobis(data.values)
    distances = pd.Series(distances, index=data.index)
    p_values = pd.Series(p_values, index=data.index)
    # Keep regions with p-value > pcut
    valid_indices = p_values[p_values > pcut].index
    filtered_regions = []
    for index in valid_indices:
        chrom = index.split(":")[0]
        start_coord = int(index.split(":")[1].split("-")[0])
        end_coord = int(index.split(":")[1].split("-")[1])
        filtered_regions.append([chrom, start_coord, end_coord])
    return filtered_regions


######################################
# Normalization Functions            #
######################################
def _norm(bw_filepath, chrom, noise, bg_sf, fg_sf, sig_sf_params, fout):
    """
    Normalize signal values for a given chromosome and write results to an output bedGraph file.

    :param bw_filepath: str, path to bigWig file.
    :param chrom: str, chromosome name.
    :param gmm: trained Gaussian Mixture Model.
    :param gmm_class_mapping: dict mapping GMM predicted class to [noise, signal].
    :param bg_sf: float, background scaling factor.
    :param fg_sf: float, reference scaling factor.
    :param sig_sf_params: list [alpha, beta] for signal scaling.
    :param fout: str, output file path for bedGraph.
    """
    bwi = pyBigWig.open(bw_filepath)
    intervals = bwi.intervals(chrom)
    with open(fout, "w") as fo:
        for interval in intervals:
            v = interval[-1]
            if v == 0.0:
                continue
            if v <= noise:  # noise
                v = v * bg_sf
            else:  # signal
                v = 2**(np.log2(v * fg_sf) * sig_sf_params[0] +
                        sig_sf_params[1])
            start = interval[0]
            end = interval[1]
            fo.write(f"{chrom}\t{start}\t{end}\t{v}\n")
    bwi.close()


def normTgtBw(bw_filepath,
              noise,
              bg_sf,
              fg_sf,
              sig_sf_params,
              fnOut,
              n_jobs=2,
              csf=""):
    """
    Normalize the target sample bigWig file and convert the output bedGraph to bigWig format.

    :param bw_filepath: str, target sample bigWig file.
    :param bg_sf: float, background scaling factor.
    :param sig_sf_params: list [alpha, beta] for signal scaling.
    :param fnOut: str, output file prefix.
    :param n_jobs: int, number of CPUs for parallel processing.
    :param csf: str, chromosome size file (used for bigWig conversion).
    """
    temp_dir = os.path.join(os.path.dirname(fnOut), str(random.random()))
    os.makedirs(temp_dir, exist_ok=True)

    bwi = pyBigWig.open(bw_filepath)
    chroms = list(bwi.chroms().keys())
    bwi.close()

    Parallel(n_jobs=n_jobs, backend="multiprocessing")(
        delayed(_norm)(bw_filepath, chrom, noise, bg_sf, fg_sf, sig_sf_params,
                       os.path.join(temp_dir, f"{chrom}.bdg"))
        for chrom in tqdm(chroms))
    bdg_files = sorted(glob(os.path.join(temp_dir, "*.bdg")))
    os.system("cat %s > %s.bdg" % (" ".join(bdg_files), fnOut))
    os.system(f"bedGraphToBigWig {fnOut}.bdg {csf} {fnOut}.bw")
    os.system(f"rm -fr {temp_dir}")


######################################
# Main Workflow (CLI)                #
######################################
@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.option(
    "-r",
    required=True,
    type=str,
    help=
    "Path to a BED format file specifying the reference regions. These regions are assumed to have no change between samples. Each region must have at least three columns: chromosome, start, and end. Regions with a width less than 100 bp will be ignored. Example: ref_regions.bed"
)
@click.option(
    "-o",
    required=True,
    type=str,
    help=
    "Output file prefix for the results. The program will generate several output files using this prefix, including QC plots and normalized bigWig files. Example: results/test"
)
@click.option(
    "-c",
    required=True,
    type=str,
    help=
    "Path to the control sample bigWig file. This file should be pre-normalized to RPM (Reads Per Million) or similar normalization with total reads. Example: control_sample.bw"
)
@click.option(
    "-lc",
    default="control",
    type=str,
    help=
    "Label for the control sample. This label is used in plots and output messages. Default: 'control'"
)
@click.option(
    "-t",
    required=True,
    type=str,
    help=
    "Path to the treatment sample bigWig file. This file should be pre-normalized to RPM (Reads Per Million). Example: treatment_sample.bw"
)
@click.option(
    "-lt",
    default="trt",
    type=str,
    help=
    "Label for the treatment sample. This label is used in plots and output messages. Default: 'trt'"
)
@click.option(
    "-ext",
    default=10000,
    type=int,
    help=
    "Extension size (in base pairs) to define the region around reference centers used for classification with the Gaussian Mixture Model (GMM). For narrow peaks, a typical value is 10,000 bp. For broad peaks (e.g., H3K27me3), consider increasing this value (e.g., 50,000 bp)."
)
@click.option(
    "-mode",
    required=False,
    default="norm",
    help=
    "Specifies the method used to draw the distribution of target sample to reference sample. Available options are norm (z-score like with log2 signal) and lr (linear fitting with log2 signal). Default is norm. Chose the one that final normalized aggregated signal shows the same level.",
    type=click.Choice(["norm", "lr"]),
)
@click.option(
    "-pred",
    type=click.Tuple([float, float, float, float]),
    default=(0, 0, 0, 0),
    help=
    "Skip the modeling step and use previously estimated background/signal scaling factors (from Step 4/8) and log2 data linear fitting parameters (alpha and beta, from Step 5/8), such as those derived from spike-in data."
)
@click.option(
    "-csf",
    required=True,
    type=str,
    help=
    "Path to the chromosome size file, which is required to convert bedGraph files to bigWig format. This file can be generated using the 'fetchChromSizes' command.  Example: chrom.sizes"
)
@click.option(
    "-p",
    default=2,
    type=int,
    help=
    "Number of CPUs to be used for parallel processing. Increasing this value can reduce processing time if multiple cores are available. Default: 2."
)
def paw(r, c, t, o, lc, lt, ext, mode, pred, csf, p):
    """
    PAW: Cross-sample Epigenome Data Normalization with Internal Reference Algorithm.
    
    This script normalizes a target sample (treatment/replicate 2) to a reference sample (control/replicate 1) at base-pair resolution.

    It relies on the assumptions that:

      1. Background noise levels are similar between the samples.

      2. Specific regions (e.g., conserved CTCF sites or transcription start sites) have similar signal-to-noise ratios.
     

    Prior to running PAW, ensure that the input bigWig files are pre-normalized to Reads Per Million (RPM).
   

    Examples:

      1. Typical pair-wise comparison:

         $ paw.py -r ref_regions.bed -c control_sample.bw -t treatment_sample.bw -csf mm10.chrom.sizes -o results/test      
        

      2. Replicate alignment:

         $ paw.py -r peaks.bed -c rep1.bw -t rep2.bw -o results/test -csf mm10.chrom.sizes


      3. Estimate the fitting through spike-in data then apply      
        
         $ paw.py -r si_peaks.bed -c si_wt.bw -t si_ko.bw -o results/si -csf mm10.chrom.sizes

         #read the background scaling factor from the line of Step 4/8, say 0.734

         #read the signal region fitting parameters alpha and beta from the line of Step 5/8 , say 0.934 * log2(KO) -1.433

         $ paw.py -r peaks.bed -c wt.bw -t ko.bw -o results/data -pred 0.734 0.934 -1.433 -csf mm10.chrom.sizes

    """
    start_time = datetime.now()
    script_name = os.path.basename(__file__)
    rprint(
        f"{script_name} -r {r} -o {o} -c {c} -t {t} -lc {lc} -lt {lt} -mode {mode} -pred {pred} -ext {ext} -p {p} -csf {csf}"
    )

    # Step 0: Check input files exist
    for filepath in [r, c, t, csf]:
        if not os.path.isfile(filepath):
            rprint(f"ERROR! Input file {filepath} does not exist. Return!")
    output_dir = os.path.dirname(o)
    if output_dir != "" and not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)

    # Step 1: Read reference peaks from BED file
    refPeaks = readBed(r)

    # Step 2: Remove outliers from reference peaks using the Mahalanobis Distance test
    rprint(f"[{o}] Step 1/8: remove potential outliers with MD test")
    filtered_peaks = removeOutliers(refPeaks, c, t, pcut=0.1)

    # Step 3: Generate foreground (signal) and background regions
    fg_regions, bg_regions = getFgBgs(filtered_peaks)
    rprint(
        f"[{o}] Step 2/8: reference peaks: {len(fg_regions)}; background regions: {len(bg_regions)}"
    )

    # Step 4: Visualize original signals around reference centers
    rprint(f"[{o}] Step 3/8: check original signals")
    showSig(fg_regions,
            c,
            t,
            o + "_1_orig",
            title="original signal",
            refLabel=lc,
            tgtLabel=lt,
            ext=ext)

    # Step 5: Perform initial QC for noise and signal-to-noise ratio
    rprint(
        f"[{o}] Step 4/8: initial QC for background noise level and signal-to-noise ratio."
    )
    fgRef, fgTgt, bgRef, bgTgt = getQc(fg_regions,
                                       bg_regions,
                                       c,
                                       t,
                                       o,
                                       refLabel=lc,
                                       tgtLabel=lt)

    if pred == (0, 0, 0, 0):
        bg_scaling_factor = bgRef.mean() / bgTgt.mean()
        fg_scaling_factor = fgRef.mean() / fgTgt.mean()
        rprint(
            f"[{o}] Step 4/8: estimated background scaling factor: {bg_scaling_factor:.3f}."
        )
        rprint(
            f"[{o}] Step 4/8: estimated reference scaling factor: {fg_scaling_factor:.3f}."
        )
        # Step 6: Estimate linear fit or scaling parameters for signal regions
        rprint(f"[{o}] Step 5/8: estimate signal region fitting parameters")
        alpha, beta = estFit(fg_regions,
                             c,
                             t,
                             o,
                             fg_scaling_factor,
                             lc,
                             lt,
                             mode=mode,
                             cpu=p)
        if beta > 0:
            rprint(
                f"[{o}] Step 5/8: estimated linear fitting: log2({lt}) = {alpha:.3f} * log2({fg_scaling_factor:.3f}{lt}) + {beta:.3f}"
            )
        else:
            rprint(
                f"[{o}] Step 5/8: estimated linear fitting: log2({lt}) = {alpha:.3f} * log2({fg_scaling_factor:.3f}{lt}) {beta:.3f}"
            )
    else:
        rprint(
            f"[{o}] Step 4&5/8: use background scaling factor and signal region fitting parameters through -pred {pred}"
        )
        bg_scaling_factor, fg_scaling_factor, alpha, beta = pred
        #just show the qc plot
        estFit(fg_regions,
               c,
               t,
               o,
               fg_scaling_factor,
               lc,
               lt,
               alpha=alpha,
               beta=beta,
               mode=mode)

    # Step 7: Train GMM for classification of target sample regions
    noise = getNoiseCut(fg_regions, bg_regions, t, o, ext=ext)
    rprint(
        f"[{o}] Step 6/8: cutoff {noise} is used to classification signal vs. background."
    )

    # Step 8: Normalize the treatment sample bigWig file
    rprint(f"[{o}] Step 7/8: normalize target sample bigWig file.")
    normTgtBw(
        t,
        noise,
        bg_scaling_factor,
        fg_scaling_factor,
        [alpha, beta],
        o + "_" + lt,
        n_jobs=p,
        csf=csf)

    # Step 9: Visualize corrected signals around reference centers
    rprint(f"[{o}] Step 8/8: check normalized signals")
    showSig(fg_regions,
            c,
            o + "_" + lt + ".bw",
            o + "_5_corr",
            title="normalized signal",
            refLabel=lc,
            tgtLabel=lt,
            ext=ext)

    end_time = datetime.now()
    rprint(f"{script_name} job finished. Used time: {end_time - start_time}")


if __name__ == "__main__":
    paw()
