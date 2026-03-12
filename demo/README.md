## Ryder: Epigenome Normalization and Variable Feature Identification

Ryder is a flexible Python package for the normalization and differential analysis of epigenomic data. It leverages stable internal reference regions to correct for technical artifacts genome-wide, supporting a wide range of assays including ChIP-seq, CUT&RUN, ATAC-seq, DNase-seq, and MNase-seq.      

This document provides instructions for using the demo datasets provided with this repository to reproduce the key results presented in our manuscript.

---
### Demo Datasets Usage

The following section details the commands to run Ryder on each of the provided datasets, corresponding to the panels in Figure 1 of our manuscript.
First, download and decompress the desired demo data package in this ryder/demo folder. For example:

``` Bash
wget https://hpc.nih.gov/~caoy7/pub/9.ryder/1.Mice_DN3_GATA3_KO_DNase-seq.tar.gz 
tar -xzvf 1.Mice_DN3_GATA3_KO_DNase-seq.tar.gz
```

Then, navigate into the created directory and run the run.sh script to reporduce the results.    
Details in run.sh are following. 

``` Bash 
paw.py -r ./CTCFnoGATA3.bed -c DNase-seq_DN3_GATA3_WT.bw -t DNase-seq_DN3_GATA3_KO.bw -lc WT -lt KO -csf ../../../refData/mm10.chrom.sizes -p 10 -o normalized
patrol.py -c ./DNase-seq_DN3_GATA3_WT.bw -t ./normalized_KO.bw -lc WT -lt KO -r ./DN3_GATA3_DHS.bed -o normalized_GATA3KO -xlim -4 2 -ylim -3 3 
```

---
### Details of Data   

#### 1. DNase-seq in GATA3 KO Mouse Thymocytes

This dataset contains DNase-seq data from wild-type (WT) and GATA3 knockout (KO) mouse DN3 thymocytes. The goal is to normalize the data and identify differential accessible sites.    

Dataset: [1.Mice_DN3_GATA3_KO_DNase-seq.tar.gz](https://hpc.nih.gov/~caoy7/pub/9.ryder/1.Mice_DN3_GATA3_KO_DNase-seq.tar.gz)

#### 2. DNase-seq in BRG1-AID Mouse Fibroblasts with Spike-in

This dataset contains DNase-seq data from WT and BRG1-AID mouse fibroblasts, with human 293T cells used as a spike-in control. This example shows how to perform internal reference normalization even when spike-ins are present. 

Dataset: [2.Mice_FibroBlast_BRG1_AID_DNase-seq.tar.gz](https://hpc.nih.gov/~caoy7/pub/9.ryder/2.Mice_FibroBlast_BRG1_AID_DNase-seq.tar.gz)     

#### 3. CUT&RUN and ATAC-seq in BRG1-dTAG mESCs
This dataset contains CUT&RUN and ATAC-seq data from mouse embryonic stem cells (mESCs) where BRG1 was progressively degraded using the dTAG system 

Dataset: [3.Mice_mESC_BRG1_dTAG_CUTandRUN_ATAC-seq.tar.gz](https://hpc.nih.gov/~caoy7/pub/9.ryder/3.Mice_mESC_BRG1_dTAG_CUTandRUN_ATAC-seq.tar.gz)    
Reference: [Hagihara & Zhang et al.](https://www.nature.com/articles/s41588-025-02305-z).

#### 4. ATAC-seq in BRG1 Inhibitor-Treated Human MV411 Cells    
This dataset contains ATAC-seq data from human MV411 cells treated with a BRG1 inhibitor. It serves as another example of correcting for global changes at enhancers.  
Dataset: [4.human_MV411_BRG1_Inhibitor_ATAC-seq.tar.gz](https://hpc.nih.gov/~caoy7/pub/9.ryder/4.human_MV411_BRG1_Inhibitor_ATAC-seq.tar.gz)     
Reference: [Martin & Ablondiet al.](https://www.cell.com/cell/fulltext/S0092-8674(23)01088-7)    

#### 5. MNase-seq in BRG1 Inhibitor-Treated mESCs
This dataset contains MNase-seq data to profile nucleosome occupancy changes in mESCs after treatment with a BRG1 inhibitor. This example demonstrates background-only scaling for subtle signals.     
Dataset: [5.Mice_mESC_BRG1_Inhibitor_MNase-seq.tar.gz](https://hpc.nih.gov/~caoy7/pub/9.ryder/5.Mice_mESC_BRG1_Inhibitor_MNase-seq.tar.gz)   
Reference: [Iurlaro et al.](https://www.nature.com/articles/s41588-020-00768-w)    

#### 6. H3K27me3 ChIP-seq in EZH2 Inhibitor-Treated PC9 Cells
This dataset contains H3K27me3 and H3K4me3 ChIP-seq data from human PC9 cells treated with an EZH2 inhibitor, which causes a global decrease in H3K27me3 but not H3K4me3. Ryder correctly quantifies this change.    
Dataset: [6.human_PC90_EZH2_Inhibitor_ChIP-seq.tar.gz](https://hpc.nih.gov/~caoy7/pub/9.ryder/6.human_PC90_EZH2_Inhibitor_ChIP-seq.tar.gz)    
Reference: [Egan et al.](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0166438)     


#### 7. H3K9ac ChIP-seq in HDAC Inhibitor-Treated HeLa Cells     
This dataset contains H3K9ac ChIP-seq data from HeLa cells treated with an HDAC inhibitor (TSA), which causes a global increase in histone acetylation. The 2.5-fold spike-in condition was omitted from this demo, as its low signal and extreme ratio make it an unrealistic use case.    
Dataset: [7.human_HeLa_H3K9ac_TSA_ChIP-seq.tar.gz](https://hpc.nih.gov/~caoy7/pub/9.ryder/7.human_HeLa_H3K9ac_TSA_ChIP-seq.tar.gz)      
Reference: [Patel et al.](https://www.nature.com/articles/s41587-024-02377-y)  
