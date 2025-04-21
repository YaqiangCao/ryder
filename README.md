# Ryder: Epigenome Normalization with Internal Reference and Variable Feature Analysis

![Rescuing Your Data with Paw Patrol!](https://i0.wp.com/cardiffmummysays.com/wp-content/uploads/2016/04/Paw-Patrol-5.jpg?resize=930%2C360)   

---
## Required python packages and other tools 
All requried other packages can be installed through conda or pip. 

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
## Install 
```
git clone --depth=1 https://github.com/YaqiangCao/ryder
cd ryder
pip install ./
```

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
       -o results/test

    2. Replicate alignment:

       $ paw.py -r peaks.bed -c rep1.bw -t rep2.bw -o results/test

Options:
  -r TEXT       Path to a BED format file specifying the reference regions.
                These regions are assumed to have no change between samples.
                Each region must have at least three columns: chromosome,
                start, and end. Regions with a width less than 100 bp will be
                ignored. Example: ref_regions.bed  [required]
  -o TEXT       Output file prefix for the results. The program will generate
                several output files using this prefix, including QC plots and
                normalized bigWig files. Example: results/test  [required]
  -c TEXT       Path to the control sample bigWig file. This file should be
                pre-normalized to RPM (Reads Per Million) or similar
                normalization with total reads. Example: control_sample.bw
                [required]
  -lc TEXT      Label for the control sample. This label is used in plots and
                output messages. Default: 'control'
  -t TEXT       Path to the treatment sample bigWig file. This file should be
                pre-normalized to RPM (Reads Per Million). Example:
                treatment_sample.bw  [required]
  -lt TEXT      Label for the treatment sample. This label is used in plots
                and output messages. Default: 'trt'
  -ext INTEGER  Extension size (in base pairs) to define the region around
                reference centers used for classification with the Gaussian
                Mixture Model (GMM). For narrow peaks, a typical value is
                10,000 bp. For broad peaks (e.g., H3K27me3), consider
                increasing this value (e.g., 50,000 bp).
  -csf TEXT     Path to the chromosome size file, which is required to convert
                bedGraph files to bigWig format. This file can be generated
                using the 'fetchChromSizes' command.  Example: chrom.sizes
                [required]
  -p INTEGER    Number of CPUs to be used for parallel processing. Increasing
                this value can reduce processing time if multiple cores are
                available. Default: 2.
  -h, --help    Show this message and exit.

```


### Detect variable features 
```
Usage: patrol.py [OPTIONS]

  PATROL: Epigenome Variable Features Detection Algorithm.

  This script detects highly variable genomic features by performing either a
  two-pass Mahalanobis Distance test ('MD' mode) or Fold Change/Poisson test
  ('FC' mode) on normalized ChIP-seq/DNase-seq/ATAC-seq data by paw.py. It
  quantifies signals over specified regions, estimates background noise, and
  identifies differential signals using statistical tests. It Outputs
  statistics, MA plots, aggregate plots, and BED files of significant regions.

  Examples:

    $ patrol.py -r regions.bed -c control.bw -t treatment.bw -o results_diff

Options:
  -r TEXT        Path to a BED format file specifying regions of interest for
                 differential analysis. Each region must have at least three
                 columns: chromosome, start, and end. Example: regions.bed
                 [required]
  -o TEXT        Output file prefix. Several output files (statistics and
                 plots) will be generated using this prefix. Example:
                 results/test  [required]
  -c TEXT        Path to the control sample bigWig file. This file should be
                 normalized (e.g., to RPM). Example: control.bw  [required]
  -lc TEXT       Label for the control sample (default: 'control').
  -t TEXT        Path to the treatment sample bigWig file. This file should be
                 normalized (e.g., to RPM). Example: treatment.bw  [required]
  -lt TEXT       Label for the treatment sample (default: 'trt').
  -pcut FLOAT    P-value cutoff for the two-pass Mahalanobis distance test if
                 -mode is MD or Poisson test if -mode is FC. Regions with
                 p-values below this cutoff are considered significant
                 variable features. Default is 0.01.
  -mode [MD|FC]  Variable feature detection mode: 'MD' (two-pass Mahalanobis
                 distance test) or 'FC' (Fold Change > 2 with Poisson test).
                 Default is FC.
  -h, --help     Show this message and exit.
```
