# Ryder: Epigenome normalization using a two-tier model and internal reference regions <img src="./pawPatrol.png" width="5%" />

---
## Introduction   
Ryder is a flexible Python package for the normalization and differential analysis of epigenomic data. It leverages stable internal reference regions to correct for technical artifacts genome-wide, supporting a wide range of assays including ChIP-seq, CUT&RUN, ATAC-seq, DNase-seq, and MNase-seq.

If you find Ryder useful, please give us a star at github and cite our paper:
Yaqiang Cao, Guangzhe Ge, and Keji Zhao (2026) "Ryder: Epigenome normalization using a two-tier model and internal reference regions"

---
### Required python packages and other tools 
All required other packages can be installed through conda or pip. 

python packages:
```
numpy 
scipy
pandas
matplotlib
seaborn 
click 
joblib
tqdm
pyBigWig
scikit-learn
```

Linux command line tools:
```
bedGraphToBigWig
```
---
### Install from scratch with conda or [mamba](https://github.com/mamba-org/mamba) 
```
mamba env create -n ryder python=3.12
conda activate ryder 
mamba install numpy scipy pandas matplotlib seaborn click joblib tqdm scikit-learn 
mamba install pybigwig ucsc-bedgraphtobigwig -c bioconda -c conda-forge
git clone --depth=1 https://github.com/YaqiangCao/ryder
cd ryder
pip install ./
```

---
### Install with pre-installed packages
```
git clone --depth=1 https://github.com/YaqiangCao/ryder
cd ryder
pip install ./
```

---
## Demo      
All datasets and usage examples required to reproduce the key analyses presented in the manuscript are available in the [demo](https://github.com/YaqiangCao/ryder/tree/main/demo) directory.    

---
## Usage

### Normalize 
```
Usage: paw.py [OPTIONS]

  PAW: Cross-sample Epigenome Data Normalization with Internal Reference
  Algorithm.

  This script normalizes a target sample (treatment/replicate 2) to a
  reference sample (control/replicate 1) at base-pair resolution.

  It relies on the assumptions that:

    1. Background noise levels are similar between the samples.

    2. Specific regions (e.g., conserved CTCF sites or transcription start
    sites) have similar signal-to-noise ratios.

  Prior to running PAW, ensure that the input bigWig files are pre-normalized
  to Reads Per Million (RPM).

  Examples:

    1. Typical pair-wise comparison:

       $ paw.py -r ref_regions.bed -c control_sample.bw -t treatment_sample.bw
       -csf mm10.chrom.sizes -o results/test

    2. Replicate alignment:

       $ paw.py -r peaks.bed -c rep1.bw -t rep2.bw -o results/test -csf
       mm10.chrom.sizes

    3. Estimate the fitting through spike-in data then apply

       $ paw.py -r si_peaks.bed -c si_wt.bw -t si_ko.bw -o results/si -csf
       mm10.chrom.sizes

       #read the background scaling factor from the line of Step 4/8, say
       0.734

       #read the reference scaling factor from the line of Step 4/8, say 1.223

       #read the signal region fitting parameters alpha and beta from the line
       of Step 5/8 , say 0.934 * log2(1.223*KO) -1.433

       $ paw.py -r peaks.bed -c wt.bw -t ko.bw -o results/data -pred 0.734
       1.223 0.934 -1.433 -csf mm10.chrom.sizes

    4. only apply the scaling factor estimated from background (or peaks) such
    as the MNase-seq data

       $ paw.py -r Nuc.bed -c wt.bw -t ko.bw -o results/data -flat bg -csf
       mm10.chrom.sizes

Options:
  -r TEXT                         Path to a BED format file specifying the
                                  reference regions. These regions are assumed
                                  to have no change between samples. Each
                                  region must have at least three columns:
                                  chromosome, start, and end. Regions with a
                                  width less than 100 bp will be ignored.
                                  Example: ref_regions.bed  [required]
  -o TEXT                         Output file prefix for the results. The
                                  program will generate several output files
                                  using this prefix, including QC plots and
                                  normalized bigWig files. Example:
                                  results/test  [required]
  -c TEXT                         Path to the control sample bigWig file. This
                                  file should be pre-normalized to RPM (Reads
                                  Per Million) or similar normalization with
                                  total reads. Example: control_sample.bw
                                  [required]
  -lc TEXT                        Label for the control sample. This label is
                                  used in plots and output messages. Default:
                                  'control'
  -t TEXT                         Path to the treatment sample bigWig file.
                                  This file should be pre-normalized to RPM
                                  (Reads Per Million). Example:
                                  treatment_sample.bw  [required]
  -lt TEXT                        Label for the treatment sample. This label
                                  is used in plots and output messages.
                                  Default: 'trt'
  -ext INTEGER                    Extension size (in base pairs) to define the
                                  region around reference centers used for
                                  classification with the Gaussian Mixture
                                  Model (GMM). For narrow peaks, a typical
                                  value is 10,000 bp. For broad peaks (e.g.,
                                  H3K27me3), consider increasing this value
                                  (e.g., 50,000 bp). For some or specific
                                  DNase-seq/ATAC-seq footprinting analysis,
                                  narrow down to 200 if default one does not
                                  work well.
  -mode [norm|lr]                 Specifies the method used to draw the
                                  distribution of target sample to reference
                                  sample. Available options are norm (z-score
                                  like with log2 signal) and lr (linear
                                  fitting with log2 signal). Default is norm.
                                  Chose the one that final normalized
                                  aggregated signal shows the same level.
  -pred <FLOAT FLOAT FLOAT FLOAT>...
                                  Skip the modeling step and use previously
                                  estimated background/signal scaling factors
                                  (from Step 4/8) and log2 data linear fitting
                                  parameters (alpha and beta, from Step 5/8),
                                  such as those derived from spike-in data.
  -noise FLOAT                    Noise cutoff to override default estimating.
  -csf TEXT                       Path to the chromosome size file, which is
                                  required to convert bedGraph files to bigWig
                                  format. This file can be generated using the
                                  'fetchChromSizes' command.  Example:
                                  chrom.sizes  [required]
  -p INTEGER                      Number of CPUs to be used for parallel
                                  processing. Increasing this value can reduce
                                  processing time if multiple cores are
                                  available. Default: 2.
  -flat [none|fg|bg]              Set this flag if data is MNase-seq or in
                                  4_NoiseCutoff.pdf can not distinguish
                                  foreground signal and background noise. Fg
                                  mode will ignore the scaling factor for the
                                  background region and bg mode will only use
                                  the scaling factor estimated from background
                                  region.
  -pcut FLOAT                     P-value cutoff for the Mahalanobis distance
                                  test to remove noise. Regions with p-values
                                  below this cutoff are considered outliers
                                  and removed from scaling factor modeling.
                                  Default is 0.1.
  -h, --help                      Show this message and exit.

```

-------
-------
### Detect variable features 
```
Usage: patrol.py [OPTIONS]

  PATROL: Epigenome Variable Features Detection Algorithm.

  This script detects highly variable genomic features by performing either a
  two-pass Mahalanobis Distance test ('MD' mode), or Fold Change/Poisson test
  ('FC' mode), or just Fold Change (FCn) on normalized ChIP-seq/DNase-
  seq/ATAC-seq data by paw.py. It quantifies signals over specified regions,
  estimates background noise, and identifies differential signals using
  statistical tests. It Outputs statistics, MA plots, aggregate plots, and BED
  files of significant regions.

  Examples:

    $ patrol.py -r regions.bed -c control.bw -t treatment.bw -o results_diff

Options:
  -r TEXT                 Path to a BED format file specifying regions of
                          interest for differential analysis. Each region must
                          have at least three columns: chromosome, start, and
                          end. Example: regions.bed  [required]
  -o TEXT                 Output file prefix. Several output files (statistics
                          and plots) will be generated using this prefix.
                          Example: results/test  [required]
  -c TEXT                 Path to the control sample bigWig file. This file
                          should be normalized (e.g., to RPM). Example:
                          control.bw  [required]
  -lc TEXT                Label for the control sample (default: 'control').
  -t TEXT                 Path to the treatment sample bigWig file. This file
                          should be normalized (e.g., to RPM). Example:
                          treatment.bw  [required]
  -lt TEXT                Label for the treatment sample (default: 'trt').
  -pcut FLOAT             P-value cutoff for the two-pass Mahalanobis distance
                          test if -mode is MD or Poisson test if -mode is FC.
                          Regions with p-values below this cutoff are
                          considered significant variable features. Default is
                          0.01.
  -mode [MD|FC|FCn]       Variable feature detection mode: 'MD' (two-pass
                          Mahalanobis distance test), 'FC' (Fold Change > 2
                          with Poisson test), 'FCn' (Fold Change > 2 without
                          Poisson test). Default is FCn.
  -xlim <FLOAT FLOAT>...  Set xlim for the MA plot.
  -ylim <FLOAT FLOAT>...  Set ylim for the MA plot.
  -h, --help              Show this message and exit.

```

---
*The package's name, along with the `paw.py` and `patrol.py` scripts, was inspired by the many hours my son Will and I spent watching PAW Patrol together during its development. Whenever a job's too big, no pup's too small!*
