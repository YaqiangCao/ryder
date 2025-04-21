#python paw.py -r ../1.testGuagnzheDNaseGATA3KO/1.raw/DN3_Ref.bed -c ../1.testGuagnzheDNaseGATA3KO/1.raw/DNase-seq_DN3_GATA3_WT.bw -t ../1.testGuagnzheDNaseGATA3KO/1.raw/DNase-seq_DN3_GATA3_KO.bw -o test/test -p 10 -csf ~/caoy7/Projects/0.Reference/2.mm10/1.fa/mm10.chrom.sizes.main
#python paw.py -r ../2.testGuangzheDNaseBRG1KO/1.raw/CTCF.bed -t ../2.testGuangzheDNaseBRG1KO/1.raw/DP_BRG1_KO.bw -c ../2.testGuangzheDNaseBRG1KO/1.raw/DP_BRG1_WT.bw -o test/BRG1 -p 20 -csf ~/caoy7/Projects/0.Reference/2.mm10/1.fa/mm10.chrom.sizes.main -lt KO -lc WT 

#python patrol.py -r ../2.testGuangzheDNaseBRG1KO/1.raw/DP_BRG1_peaks.bed -c ../2.testGuangzheDNaseBRG1KO/1.raw/DP_BRG1_WT.bw -t test/BRG1_KO.bw -o test/BRG1_BRG1 -lc WT -lt BRG1KO -pcut 0.01
#python patrol.py -r ../2.testGuangzheDNaseBRG1KO/1.raw/CTCF.bed -c ../2.testGuangzheDNaseBRG1KO/1.raw/DP_BRG1_WT.bw -t test/BRG1_KO.bw -o test/BRG1_CTCF -lc WT -lt BRG1KO -pcut 0.01

#revised paw.py
#python paw.py -r ../2.testGuangzheDNaseBRG1KO/1.raw/CTCF.bed -t ../2.testGuangzheDNaseBRG1KO/1.raw/DP_BRG1_KO.bw -c ../2.testGuangzheDNaseBRG1KO/1.raw/DP_BRG1_WT.bw -o test3/BRG1 -p 20 -csf ~/caoy7/Projects/0.Reference/2.mm10/1.fa/mm10.chrom.sizes.main -lt KO -lc WT 
python patrol.py -r ../2.testGuangzheDNaseBRG1KO/1.raw/CTCF.bed -c ../2.testGuangzheDNaseBRG1KO/1.raw/DP_BRG1_WT.bw -t test3/BRG1_KO.bw -o test4/BRG1_CTCF -lc WT -lt BRG1KO &
python patrol.py -r ../2.testGuangzheDNaseBRG1KO/1.raw/CTCF.bed -c ../2.testGuangzheDNaseBRG1KO/1.raw/DP_BRG1_WT.bw -t test3/BRG1_KO.bw -o test4/BRG1_CTCF_FC -mode FC -lc WT -lt BRG1KO &
python patrol.py -r ../2.testGuangzheDNaseBRG1KO/1.raw/DP_BRG1_peaks.bed -c ../2.testGuangzheDNaseBRG1KO/1.raw/DP_BRG1_WT.bw -t test3/BRG1_KO.bw -o test4/BRG1_BRG1 -lc WT -lt BRG1KO  &
python patrol.py -r ../2.testGuangzheDNaseBRG1KO/1.raw/DP_BRG1_peaks.bed -c ../2.testGuangzheDNaseBRG1KO/1.raw/DP_BRG1_WT.bw -t test3/BRG1_KO.bw -o test4/BRG1_BRG1 -lc WT -lt BRG1KO -mode FC -pcut 0.01 &
