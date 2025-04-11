#!/usr/bin/env python
#--coding:utf-8--
"""
PAW algorithm for cross-sample ChIP-seq/DNase-seq/ATAC-seq normalization with internal control. 

2025-04-04: init
2025-04-09: revised and basically finished
"""

__author__ = "CAO Yaqiang"
__date__ = "2025-04-04"
__modified__ = ""
__email__ = "caoyaqiang0410@gmail.com"

#sys library
import os
import random
import warnings
from glob import glob
from datetime import datetime

#3rd library
import click
import pyBigWig
import numpy as np
import pandas as pd
from tqdm import tqdm
from joblib import Parallel, delayed
from sklearn.mixture import GaussianMixture as GMM

#global settings
warnings.filterwarnings("ignore")
#plotting setting
import matplotlib as mpl

mpl.use("pdf")
import seaborn as sns

mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["figure.figsize"] = (3.2, 2.2)
mpl.rcParams["savefig.transparent"] = True
mpl.rcParams["savefig.bbox"] = "tight"
mpl.rcParams["font.size"] = 7.0
mpl.rcParams["font.sans-serif"] = "Arial"
mpl.rcParams["savefig.format"] = "pdf"
import pylab

sns.set_style("white")
colors = [
    (0.8941176470588236, 0.10196078431372549, 0.10980392156862745),
    (0.21568627450980393, 0.49411764705882355, 0.7215686274509804),
    (0.30196078431372547, 0.6862745098039216, 0.2901960784313726),
    (0.596078431372549, 0.3058823529411765, 0.6392156862745098),
    (1.0, 0.4980392156862745, 0.0),
    (0, 0, 0),  #black
    #(1.0, 1.0, 0.2),
    #(0.7, 1.0, 0.2),
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


def rprint(message):
    """
    Print message with time.
    @param message: str, 
    @returns: None
    """
    report = "\t".join([
        str(datetime.now()),
        message,
    ])
    print(report)


def readBed(f):
    rs = []
    for line in open(f):
        line = line.split("\n")[0].split("\t")
        if len(line) < 3:
            continue
        try:
            line[1] = int(line[1])
            line[2] = int(line[2])
        except:
            continue
        if line[2] - line[1] < 100:
            continue
        rs.append(line[:3])
    return rs


def showSig(rs,
            cf,
            tf,
            o,
            title="",
            refLabel="ref",
            tgtLabel="tgt",
            ext=10000):
    x = np.arange(-ext, ext + 1)
    cs = np.zeros(len(x))
    ts = np.zeros(len(x))
    cfo = pyBigWig.open(cf)
    tfo = pyBigWig.open(tf)
    for r in tqdm(rs):
        chrom = r[0]
        center = int((r[1] + r[2]) / 2)
        sa = cfo.values(chrom, center - ext, center + ext + 1)
        sa = np.nan_to_num(sa)
        cs += sa
        sb = tfo.values(chrom, center - ext, center + ext + 1)
        sb = np.nan_to_num(sb)
        ts += sb
    cs = cs / len(rs)
    ts = ts / len(rs)
    fig, ax = pylab.subplots()
    ax.plot(x, cs, label=refLabel, color=colors[0])
    ax.plot(x, ts, label=tgtLabel, color=colors[1])
    ax.set_ylabel("avg. signal")
    ax.set_xlabel("distance from reference center")
    ax.set_title(title)
    ax.set_xlim([-int(ext / 2), int(ext / 2) + 1])
    leg = ax.legend(labelcolor=[colors[0], colors[1]],
                    frameon=False,
                    markerscale=0)
    pylab.savefig(o + ".pdf")


def buildCov(rs):
    """
    Genomic region coverage
    """
    lims = {}
    cov = {}
    for r in rs:
        #coverages
        if r[0] not in cov:
            cov[r[0]] = set()
        cov[r[0]].update(range(r[1], r[2] + 1))
        #range limitations
        if r[0] not in lims:
            lims[r[0]] = [r[1], r[2]]
        if r[1] < lims[r[0]][0]:
            lims[r[0]][0] = r[1]
        if r[2] > lims[r[0]][1]:
            lims[r[0]][1] = r[2]
    return cov, lims


def getRegion(cov, margin=1):
    """
    Get regions from coverage
    """
    rs = []
    for c, s in cov.items():
        s = list(s)
        s.sort()
        i = 0
        while i < len(s) - 1:
            for j in range(i + 1, len(s)):
                if s[j] - s[j - 1] > margin:
                    break
                else:
                    continue
            start = s[i]
            end = s[j - 1]
            i = j  #update search start
            rs.append([c, start, end])
    return rs


def checkBgOverlaps(c, s, e, cov, lims):
    if s < lims[c][0] or e > lims[c][1]:
        return True
    for i in range(s, e + 1):
        if i in cov[c]:
            return True
    return False


def getFgBgs(refPeaks, exts=[5, 10, 20]):
    """
    Get the background regions.
    """
    refCov, refLims = buildCov(refPeaks)
    fgs = getRegion(refCov)  #sorted reference
    bgs = []
    for r in fgs:
        d = r[2] - r[1]
        c = r[0]
        for ext in exts:
            s = r[1] - ext * d
            e = r[2] - ext * d
            if checkBgOverlaps(c, s, e, refCov, refLims):
                continue
            bgs.append([c, s, e])
            s = r[1] + ext * d
            e = r[2] + ext * d
            if checkBgOverlaps(c, s, e, refCov, refLims):
                continue
            bgs.append([c, s, e])
    return fgs, bgs


def getBinMean(s, bins=100):
    """
    Get the mean of bins for a array.
    @param s: np.array
    @param bins: int, how many bins as converted
    """
    width = int(len(s) / bins)
    ns = s[:bins * width].reshape(-1, width).mean(axis=1)
    return ns


def _getBwSig(r, bw, bins=100):
    """
    Get the signal from bigWig file and get the bin mean.
    """
    bwo = pyBigWig.open(bw)
    try:
        ns = bwo.values(r[0], r[1], r[2])
        ns = np.nan_to_num(ns)
        if len(ns) > bins:
            ns = getBinMean(ns, bins)
            return ns
        else:
            return np.zeros(bins)
    except:  #fetching data error
        return np.zeros(bins)


def getBwSig(rs, f, bins=100, p=2):
    """
    Get region signal from bigWig file. 
    """
    s = Parallel(n_jobs=p,
                 backend="multiprocessing")(delayed(_getBwSig)(r, f, bins)
                                            for r in rs)
    return np.array(s)


def plotBinSig(ax, x, y1, y2, label1, label2, xlabel, ylabel, title):
    ax.plot(x, y1, label=label1, color=colors[0])
    ax.plot(x, y2, label=label2, color=colors[1])
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    leg = ax.legend(labelcolor=[colors[0], colors[1]],
                    frameon=False,
                    markerscale=0)
    return ax


def getQc(fgs,
          bgs,
          refBw,
          tgtBw,
          fnOut,
          refLabel="ref",
          tgtLabel="tgt",
          bins=100):
    """
    Quality control: 
    1. noise compare
    2. signal to noise ratio compare
    """
    fgRef = getBwSig(fgs, refBw, bins=bins)
    fgRef = fgRef.mean(axis=0)
    fgTgt = getBwSig(fgs, tgtBw, bins=bins)
    fgTgt = fgTgt.mean(axis=0)
    bgRef = getBwSig(bgs, refBw, bins=bins)
    bgRef = bgRef.mean(axis=0)
    bgTgt = getBwSig(bgs, tgtBw, bins=bins)
    bgTgt = bgTgt.mean(axis=0)
    bgfc = bgRef.mean() / bgTgt.mean()
    fgfc = fgRef.mean() / fgTgt.mean()

    #show the results
    fig, axs = pylab.subplots(1, 3, figsize=(6, 2), sharex=True)
    axs = axs.reshape(-1)
    x = np.arange(bins)
    #background level
    ax = axs[0]
    plotBinSig(ax, x, bgRef, bgTgt, refLabel, tgtLabel, "bins", "avg. signal",
               "background\n sf:%.3f" % bgfc)
    #signal level
    ax = axs[1]
    plotBinSig(ax, x, fgRef, fgTgt, refLabel, tgtLabel, "bins", "avg. signal",
               "reference region\n sf:%.3f" % fgfc)
    #signal to noise ratio
    refSN = fgRef / bgRef
    tgtSN = fgTgt / bgTgt
    ax = axs[2]
    snfc = refSN.mean() / tgtSN.mean()
    plotBinSig(ax, x, refSN, tgtSN, refLabel, tgtLabel, "bins", "S2N ratio",
               "signal to noise ratio\n sf:%.3f" % snfc)
    pylab.tight_layout()
    pylab.savefig(fnOut + "_2_qc.pdf")
    return fgRef, fgTgt, bgRef, bgTgt


def estFit(fgs,
           refBw,
           tgtBw,
           bgRef,
           bgTgt,
           sf,
           fnOut,
           refLabel="ref",
           tgtLabel="tgt",
           bins=2):
    """
    Estimate linear fiting between samples in the signal regions
    """
    refS = getBwSig(fgs, refBw, bins=bins)
    refS = pd.Series(refS.reshape(-1))
    tgtS = getBwSig(fgs, tgtBw, bins=bins)
    tgtS = pd.Series(tgtS.reshape(-1))

    fig, axs = pylab.subplots(1, 3, figsize=(6, 2))
    axs = axs.reshape(-1)

    #signal conversion
    """
    refS = refS - bgRef.mean()
    tgtS = (tgtS - bgTgt.mean()) * sf
    """
    s = refS[refS > 0].index
    s = tgtS[s]
    s = s[s > 0].index
    refS = refS[s]
    tgtS = tgtS[s]
    #log transformation
    refS = np.log2(refS)
    tgtS = np.log2(tgtS)

    #plot signal distribution
    ax = axs[0]
    sns.kdeplot(refS, label=refLabel, ax=ax, fill=True, color=colors[0])
    sns.kdeplot(tgtS, label=tgtLabel, ax=ax, fill=True, color=colors[1])
    ax.legend()
    ax.set_xlabel("log2(signal)")
    ax.set_title("signal distribution")

    #plot signal linear fitting
    ax = axs[1]
    ax.scatter(refS, tgtS, s=0.1)
    _m = tgtS - refS
    _a = (refS + tgtS) / 2
    _pcc = _m.corr(_a)
    s = np.min([np.min(refS), np.min(tgtS)])
    e = np.max([np.max(refS), np.max(tgtS)])
    ax.plot([s, e], [s, e], color="gray", linestyle="--")
    ax.set_title(f"before correction \n M~A PCC:{_pcc:.3f}")
    ax.set_xlabel(refLabel)
    ax.set_ylabel(tgtLabel)

    #distribution match
    ax = axs[2]
    #tgtSc = (tgtS - tgtS.mean())/tgtS.std()*refS.std() + refS.mean()
    alpha = refS.std() / tgtS.std()
    beta = refS.mean() - alpha * tgtS.mean()
    tgtSc = tgtS * alpha + beta

    m = tgtSc - refS
    a = (refS + tgtSc) / 2
    ax.scatter(a, m, s=0.1)
    if beta > 0:
        ax.set_title("after correction M~A PCC:%.3f\n%s=%.3f%s+%.3f" %
                     (m.corr(a), tgtLabel, alpha, tgtLabel, beta))
    else:
        ax.set_title("after correction M~A PCC:%.3f\n(%s)=%.3f(%s)%.3f" %
                     (m.corr(a), tgtLabel, alpha, tgtLabel, beta))

    ax.axhline(0, color="gray", linestyle="--")
    ax.set_xlabel("A, (log2(%s)+log2(%s))/2)" % (refLabel, tgtLabel))
    ax.set_ylabel("M, log2(%s)-log2(%s)" % (tgtLabel, refLabel))

    pylab.tight_layout()
    pylab.savefig(fnOut + "_3_fgSignalConversion.pdf")
    return [alpha, beta]



def getGmm(fg, bg, bw, ax, title, ext=5000, bins=10):
    """
    Train one GMM. Data should log2 first
    """
    rs = []
    for r in fg:
        center = int((r[1] + r[2]) / 2)
        nr = [r[0], center - ext, center + ext]
        rs.append(nr)
    #foreground signal
    fgS = getBwSig(fg, bw, bins=bins)
    fgS = fgS.reshape(-1)
    fgS = np.log2(fgS[fgS > 0])
    #background signal
    bgS = getBwSig(bg, bw, bins=bins)
    bgS = bgS.reshape(-1)
    bgS = np.log2(bgS[bgS > 0])
    #peak nearby region
    mixS = getBwSig(rs, bw, bins=bins)
    mixS = mixS.reshape(-1)
    mixS = np.log2(mixS[mixS > 0])
    #plot
    sns.kdeplot(bgS, label="background", fill=False, ax=ax)
    sns.kdeplot(fgS, label="reference region", fill=False, ax=ax)
    sns.kdeplot(mixS, label="reference region nearby", fill=False, ax=ax)
    #train gmm
    gmm = GMM(n_components=2,
              covariance_type="full",
              random_state=123,
              means_init=[[np.mean(fgS)],
                          [np.mean(bgS)]]).fit([[v] for v in mixS])
    ms = gmm.means_.reshape(-1)
    ws = gmm.weights_
    ax.legend()
    ax.set_xlabel("log2(signal)")
    ax.set_title(title + "\n" + "means:%.3f, %.3f\nweights:%.3f, %.3f" %
                 (ms[0], ms[1], ws[0], ws[1]))
    #correct gmm predict targets, 0 as noise , 1 as signal
    if ms[0] > ms[1]:
        cs = {0: 1, 1: 0}
    else:
        cs = {0: 0, 1: 1}
    return gmm, cs


def trainGmm(
    fgs,
    bgs,
    tgtBw,
    fnOut,
    tgtLabel="tgt",
    ext=5000,
):
    """
    Train GMM to classify background regions or signal regions
    """
    fig, ax = pylab.subplots()
    tgtGmm, tgtCs = getGmm(fgs, bgs, tgtBw, ax, tgtLabel, ext)
    #pylab.tight_layout()
    pylab.savefig(fnOut + "_4_GMM.pdf")
    return tgtGmm, tgtCs


def _norm(bw, chrom, gmm, gmmCs, noise, sf, sf2, fout):
    bwi = pyBigWig.open(bw)
    ts = bwi.intervals(chrom)
    with open(fout, "w") as fo:
        for s in ts:
            v = s[-1]
            if v == 0:
                continue
            #gmm was trained with log2 data
            t = gmmCs[gmm.predict([[np.log2(v)]])[0]]
            if t == 0:  # noise
                v = v * sf
            else:  #signal
                v = 2**(np.log2(v) * sf2[0] + sf2[1])
            start = s[0]
            end = s[1]
            line = f"{chrom}\t{start}\t{end}\t{v}\n"
            fo.write(line)
    bwi.close()


def normTgtBw(bw, gmm, gmmCs, noise, sf, sf2, fnOut, p=2, csf=""):
    """
    Normalize target sample bigWig file.
    """
    td = "/".join(fnOut.split("/")[:-1]) + "/" + str(random.random())
    os.makedirs(td)
    bwi = pyBigWig.open(bw)
    chroms = bwi.chroms().keys()
    bwi.close()
    Parallel(n_jobs=p, backend="multiprocessing")(
        delayed(_norm)(bw, chrom, gmm, gmmCs, noise, sf, sf2, td + "/" +
                       chrom + ".bdg") for chrom in tqdm(chroms))
    fs = glob(td + "/*.bdg")
    fs.sort()
    c1 = "cat %s > %s.bdg" % (" ".join(fs), fnOut)
    c2 = f"bedGraphToBigWig {fnOut}.bdg {csf} {fnOut}.bw"
    c3 = f"rm -fr {td}"
    for c in [c1, c2, c3]:
        os.system(c)


@click.command()
@click.option(
    "-r",
    required=True,
    help=
    "A BED format file specifying the reference regions, which are assumed to exhibit no change between samples.",
    type=str,
)
@click.option(
    "-o",
    required=True,
    help="Output file prefix",
    type=str,
)
@click.option(
    "-c",
    required=True,
    help="The bigWig file for the control sample.",
    type=str,
)
@click.option(
    "-lc",
    required=False,
    help="The label for the reference sample. Default is control.",
    type=str,
    default="control",
)
@click.option(
    "-t",
    required=True,
    help="The bigWig file for the treatment sample.",
    type=str,
)
@click.option(
    "-lt",
    required=False,
    help="The label for the treatment sample. Default is trt.",
    type=str,
    default="trt",
)
@click.option(
    "-ext",
    required=False,
    help=
    "To classify background and signal regions, a Gaussian Mixture Model is built using data extended from reference region centers. This parameter specifies the extension size in base pairs (default: 10,000 bp). For broad peaks, like H3K27me3, increase this value to capture enough surrounding regions (e.g., 50,000 bp).",
    type=int,
    default=10000,
)
@click.option(
    "-csf",
    required=True,
    help=
    "Chromosome size file. Can be obtained through the command of fetchChromSizes. If given, will convert the output file to bigWig.",
    type=str,
)
@click.option("-p",
              default=2,
              type=int,
              help="Number of CPUs to finish the job, default is set to 2.")
def paw(r, c, t, o, lc, lt, ext=10000, csf="", p=2):
    """
    To normalize target sample ChIP-seq, ATAC-seq, or DNase-seq data to a reference sample at base-pair resolution, we assume: 1) background noise levels are similar; and 2) specific regions, such as conserved CTCF sites or transcription start sites (TSS) of genes with unchanged expression, maintain similar signal-to-noise ratios. This normalization aims to account for inter-sample variability. It is crucial that input BigWig files are pre-normalized to Reads Per Million (RPM) before applying this method.
    
    Examples:

        1. typical pair-wise comparsion   

            paw.py -r ref.bed -c control.bw -t trt.bw -o test 
    

        2. replicates alignment 

            paw.py -r peaks.bed -c rep1.bw -t rep2.bw -o test 
    """
    #start
    start = datetime.now()
    script = os.path.basename(__file__)
    rprint(
        f"{script} -r {r} -o {o} -c {c} -t {t} -lc {lc} -lt {lt} -ext {ext} -p {p} -csf {csf}"
    )

    #step 0 parameters check
    for f in [r, c, t, csf]:
        if not os.path.isfile(f):
            rprint(f"ERROR! Input file {f} not exits. Return!")
    od = os.path.dirname(o)
    if not os.path.exists(od):
        os.makedirs(od, exist_ok=True)

    #step 1 read reference peaks/regions
    refPeaks = readBed(r)

    #step 2 generate background region
    fgs, bgs = getFgBgs(refPeaks)
    rprint("[%s] Step 1: reference peaks: %s; background regions: %s" %
           (o, len(fgs), len(bgs)))

    #step 3 show the orignal signal around reference centers
    rprint(f"[{o}] Step 2: check original signals")
    showSig(fgs,
            c,
            t,
            o + "_1_orig",
            title="original signal",
            refLabel=lc,
            tgtLabel=lt,
            ext=ext)

    #step 4 qc for signal to noise ratio and noise level
    rprint(
        f"[{o}] Step 3: initial QC for background noise level and signal-to-noise ratio."
    )
    fgRef, fgTgt, bgRef, bgTgt = getQc(fgs,
                                       bgs,
                                       c,
                                       t,
                                       o,
                                       refLabel=lc,
                                       tgtLabel=lt)
    #scaling factor for background region
    sf = bgRef.mean() / bgTgt.mean()
    rprint(f"[{o}] Step 3: estimated background scaling factor: {sf:.3f}.")

    #step 5  estimate sample-wise fitting parameters, scaling factor for signal regions
    rprint(f"[{o}] Step 4: estimate signal region fitting parameters ")
    alpha, beta = estFit(fgs, c, t, bgRef, bgTgt, sf, o, lc, lt)
    if beta > 0:
        rprint(
            f"[{o}] Step 4: estimated linear fitting for signal region: log2({lt})={alpha:.3f}log2({lt}) + {beta:.3f}"
        )
    else:
        rprint(
            f"[{o}] Step 4: estimated linear fitting for signal region: log2({lt})={alpha:.3f}log2({lt}) {beta:.3f}"
        )

    #step 6 train GMM with fg and bg data for classifiy fg and bg regions
    rprint(
        f"[{o}] Step 5: train Gaussian Mixture Model for classfication of background and signal regions"
    )
    tgtGmm, tgtCs = trainGmm(fgs, bgs, t, o, lt, ext=ext)

    #step 7 performe normalization to tgt bigwig files
    noiseTgt = bgTgt.mean()
    rprint(f"[{o}] Step 6: normalize target sample bigWig file.")
    normTgtBw(t,
              tgtGmm,
              tgtCs,
              noiseTgt,
              sf, [alpha, beta],
              o + "_" + lt,
              p=p,
              csf=csf)

    #step 8 show the corrected signal around reference centers
    rprint(f"[{o}] Step 7: check corrected signals")
    showSig(fgs,
            c,
            o + "_" + lt + ".bw",
            o + "_5_corr",
            title="correct signal",
            refLabel=lc,
            tgtLabel=lt,
            ext=ext)

    #finished
    end = datetime.now()
    usedTime = end - start
    rprint(f"{script} job finished. Used time: {usedTime}")


if __name__ == "__main__":
    paw()
