#!/usr/bin/env python
#--coding:utf-8--
"""
PATROL algorithm for cross-sample normalized ChIP-seq/DNase-seq/ATAC-seq differential peaks detection. 
"""

#sys
__author__ = "CAO Yaqiang"
__date__ = "2025-04-09"
__modified__ = ""
__email__ = "caoyaqiang0410@gmail.com"

#sys library
import os
import warnings
from datetime import datetime

#3rd
import click
import pyBigWig
import numpy as np
import pandas as pd
from tqdm import tqdm
from scipy import linalg
from scipy.stats import chi2
from joblib import Parallel, delayed

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


def getFgBgs(refPeaks, exts=[10]):
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


def quant(rs, bgs, bwf):
    """
    Quantify genomic features from bigWig file. 

    @param rs: list, [[chrom, start,end]], region of interest
    @param bwf: str, bigWig file path
    """
    bwo = pyBigWig.open(bwf)
    noise = []
    for r in bgs:
        chrom = r[0]
        start = r[1]
        end = r[2]
        ns = bwo.values(chrom, start, end)
        ns = np.nan_to_num(ns)
        noise.append(np.sum(ns) / len(ns))
    noise = np.mean(noise)
    s = {}
    for r in tqdm(rs):
        chrom = r[0]
        start = r[1]
        end = r[2]
        rid = f"{chrom}:{start}-{end}"
        ns = bwo.values(chrom, start, end)
        ns = np.nan_to_num(ns)
        ns = np.sum(ns) / len(ns)
        if ns == 0:
            ns = noise
        s[rid] = ns
    bwo.close()
    s = pd.Series(s)
    return s, noise


def mahalanobis(mat):
    """
    Caculate the mahalanobis distance.

    according to: https://www.statology.org/mahalanobis-distance-python/

    @param mat: np.array, row is each item and column is each variable 

    @return dis: np.array,Mahalanobis distance
    @return ps: np.array, chi-square test p-values
    """
    #covariance matrix
    cov = np.cov(mat, rowvar=False)
    #inverse covariance matrix
    invCov = np.linalg.inv(cov)
    #center
    center = np.mean(mat, axis=0)
    #mahalanobis distance
    mu = mat - center
    dis = np.dot(np.dot(mu, invCov), mu.T).diagonal()
    #Chi-square test p-values for detecting outliers
    ps = 1 - chi2.cdf(dis, mat.shape[1] - 1)
    return dis, ps


def twoPassesMDTest(data, pcut=0.01):
    """
    Perform MD test with two passes. First pass use all data, second pass using cov and center from data without outliers. 
    @param data: pd.DataFrame, row is item and column is sample.
    @param pcut: Chi-Square p-value cutoffs
    """
    #first pass test
    #mahalanobis distance and pvaues
    dis, ps = mahalanobis(data.values)
    dis = pd.Series(dis, index=data.index)
    ps = pd.Series(ps, index=data.index)
    inds = ps[ps < pcut].index

    #second pass test, only using data without outliers to caculate cov matrix and center
    ndata = data.drop(inds)
    #covariance matrix
    cov = np.cov(ndata.values, rowvar=False)
    #inverse covariance matrix
    invCov = np.linalg.inv(cov)
    #center
    center = np.mean(ndata.values, axis=0)

    #mahalanobis distance for all data
    mu = data.values - center
    dis = np.dot(np.dot(mu, invCov), mu.T).diagonal()
    #Chi-square test p-values for detecting outliers from all data
    ps = 1 - chi2.cdf(dis, data.shape[1] - 1)

    dis = pd.Series(dis, index=data.index)
    ps = pd.Series(ps, index=data.index)
    return dis, ps


def showMA( m, a,lc, lt, ps,o, noise, pcut=0.01):
    """
    Show MA plots.
    """
    fig, ax = pylab.subplots()
    ax.scatter(a, m, s=0.5, color="gray", alpha=0.6, label="total %s regions" % len(a))
    nps = ps[ps<=pcut].index
    ma = a[nps]  
    ma = ma[ma>noise].index
    ms = m[ma]
    trs = ms[ms>0].index
    crs = ms[ms<0].index
    ax.scatter(a[crs], m[crs], s=1, color=colors[0], alpha=0.8, label="%s %s regions" % (len(crs), lc ))
    ax.scatter(a[trs], m[trs], s=1, color=colors[1], alpha=0.8, label="%s %s regions" % (len(trs), lt ))
    leg = ax.legend( loc='best', labelcolor=["gray", colors[0], colors[1]])
    for h in leg.legendHandles:
        h._sizes = [10]
    ax.axhline(0, color="gray", linestyle="--")
    ax.set_title(f"signal comparsion\nMahalanobis distance P-value < {pcut}")
    ax.set_xlabel(f"log2({lt})+log2({lc})")
    ax.set_ylabel(f"log2({lt})-log2({lc})")
    pylab.savefig(o+"_MA.pdf")



@click.command()
@click.option(
    "-r",
    required=True,
    help=
    "A BED format file specifying the regions of interest to detect changes.",
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
    "-pcut",
    required=False,
    help=
    "The p-value cutoff for Chi-squared test for two-pass Mahalanobis distance. Default is 0.01.",
    type=float,
    default="0.05",
)
def patrol(r, c, t, o, lc, lt, pcut=0.05):
    """
    Identify highly variable features with two passes MD test after normalization with internal reference from DNase-seq, ATAC-seq, ChIP-seq and MNase-seq. 
    
    Examples:

        1. typical pair-wise comparsion   

            patrol.py -r regionOfInterest.bed -c control.bw -t pawed_trt.bw -o test 
    
    """
    #start
    start = datetime.now()
    script = os.path.basename(__file__)
    rprint(f"{script} -r {r} -o {o} -c {c} -t {t} -lc {lc} -lt {lt}")

    #step 0 parameters check
    for f in [r, c, t]:
        if not os.path.isfile(f):
            rprint(f"ERROR! Input file {f} not exits. Return!")
    od = os.path.dirname(o)
    if not os.path.exists(od):
        os.makedirs(od, exist_ok=True)

    #step 1 read reference peaks/regions; background regions were used to get the noise level to avoid signal region log2(0)
    rs = readBed(r)
    rs, bgs = getFgBgs(rs)

    #step 2 quantify features
    rprint(f"[{o}] Step 1: quantify {lc} sample")
    sc,cn = quant(rs, bgs, c)
    rprint(f"[{o}] Step 1: quantify {lt} sample")
    st,tn = quant(rs, bgs, t)

    #step 3 perform two-passes MD test
    rprint(f"[{o}] Step 2: perform the two-passes MD test")
    sc2 = np.log2(sc)
    st2 = np.log2(st)
    a = (sc2 + st2) / 2
    m = st2 - sc2
    ds = pd.DataFrame({"m": m, "a": a})
    dis, ps = twoPassesMDTest(ds, pcut=pcut)
    data = {lc: sc, lt:st, f"log2 fold change ({lt}/{lc})":m, "average": a,"MD":dis, "Chi-squared P-value": ps}
    data = pd.DataFrame(data)
    data.to_csv(f"{o}_stat.txt",sep="\t",index_label="regionId")
    rprint(f"[{o}] Step 2: test results saved to {o}_stat.txt.")
    
    #step 4 show the results
    rprint(f"[{o}] Step 3: show the detected variables MA plot")
    noise = np.log2(cn) + np.log2(tn)
    #_sc = sc[sc>cn]
    #_st = st[st>tn]
    #noise = (np.log2(np.mean(_sc)) + np.log2(np.mean(_st)))/2
    showMA(m, a, lc, lt, ps, o, noise, pcut=pcut)


    #finished
    end = datetime.now()
    usedTime = end - start
    rprint(f"{script} job finished. Used time: {usedTime}")


if __name__ == "__main__":
    patrol()
