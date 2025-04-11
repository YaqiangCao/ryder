#python paw.py -r ../1.testGuagnzheDNaseGATA3KO/1.raw/DN3_Ref.bed -c ../1.testGuagnzheDNaseGATA3KO/1.raw/DNase-seq_DN3_GATA3_WT.bw -t ../1.testGuagnzheDNaseGATA3KO/1.raw/DNase-seq_DN3_GATA3_KO.bw -o test/test -p 10 -csf ~/caoy7/Projects/0.Reference/2.mm10/1.fa/mm10.chrom.sizes.main
#python paw.py -r ../2.testGuangzheDNaseBRG1KO/1.raw/CTCF.bed -t ../2.testGuangzheDNaseBRG1KO/1.raw/DP_BRG1_KO.bw -c ../2.testGuangzheDNaseBRG1KO/1.raw/DP_BRG1_WT.bw -o test/BRG1 -p 20 -csf ~/caoy7/Projects/0.Reference/2.mm10/1.fa/mm10.chrom.sizes.main -lt KO -lc WT 

python patrol.py -r ../2.testGuangzheDNaseBRG1KO/1.raw/DP_BRG1_peaks.bed -c ../2.testGuangzheDNaseBRG1KO/1.raw/DP_BRG1_WT.bw -t test/BRG1_KO.bw -o test/BRG1_diff -lc WT -lt BRG1KO
