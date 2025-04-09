#!/usr/bin/env python
#--coding:utf-8--
"""
PATROL algorithm for cross-sample normalized ChIP-seq/DNase-seq/ATAC-seq differential peaks detection. 
"""

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
@click.option("-p",
              default=2,
              type=int,
              help="Number of CPUs to finish the job, default is set to 2.")
def patrol(r, c, t, o, lc, lt, ext=10000, csf="", p=2):
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
    for f in [r, c, t]:
        if not os.path.isfile(f):
            rprint(f"ERROR! Input file {f} not exits. Return!")
    od = os.path.dirname(o)
    if not os.path.exists(od):
        os.makedirs(od, exist_ok=True)

    #step 1 read reference peaks/regions
    refPeaks = readBed(r)

    #finished
    end = datetime.now()
    usedTime = end - start
    rprint(f"{script} job finished. Used time: {usedTime}")


if __name__ == "__main__":
    patrol()
