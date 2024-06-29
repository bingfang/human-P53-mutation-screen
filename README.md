---
Title: "Human_P53_screen_20240628"
Author: "Bingfang Ruth Xu"
Date: "2024-06-28"

---

## Trimmatic lines


```
# paired end trimming WITHOUT trimming start seqence
java -jar /Users/xubr/anaconda3/share/trimmomatic-0.39-2/trimmomatic.jar PE -phred33 data/HuP53-20240626_S1_L001_R1_001.fastq.gz data/HuP53-20240626_S1_L001_R2_001.fastq.gz data/P53_S1_L001_R1_trimmed.fastq P53_S1_L001_R1_unpair_trimmed.fastq data/P53_S1_L001_R2_trimmed.fastq P53_S1_L001_R2_unpair_trimmed.fastq TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36


```

## Subset the fastq files


```
# if reads per samples is more than 100000, subset the fastq file to reduce computation time.
/Users/xubr/local_projects/seqtk/seqtk sample -s100 P53_samples_R1_trimmed.fastq 100000 > P53_samples_R1_sub_trimmed.fastq
/Users/xubr/local_projects/seqtk/seqtk sample -s100 P53_samples_R2_trimmed.fastq 100000 > P53_samples_R2_sub_trimmed.fastq
```

## Debarcode, Genptyping, Caculation of mutation rates


```
# input sample barcode info, then run debarcode script.
# check if P53*trimmed.fastq are only input fastq
python3 src/_0debarcode_Miseq_20190515.py

# check the filter for mutation rate, then run genotyping script.
# Open TP53 snapgene, and align sequences from "Check_seq_all.txt"
# remove low coverage samples
python3 src/_1_humanP53_genotyped_20210730.py

# input ~10bp WT, mutated sequences, and run the script of mutation rate calculation.
python3 src/_2CheckSNP_Count1.py
```
