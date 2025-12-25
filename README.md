# ChIP-seq_Flow_Athaliana

This repository provides a streamlined workflow for analyzing ChIP-seq data, leveraging Python and R tools to handle data preprocessing, visualization, and downstream analysis. Follow the instructions below to set up the environment and get started.
1. merge 4 separated fastq files into one.
2. mapping with bowtie2 using Athaliana index file
3. binaryize and sort with samtools
4. select by macs3 peakcall with significant (p < 0.05) for treated against control
5. annotate peaks with neighboring genes using R ChIPpeakAnno

---

## Prerequisites

### Required Software
- Python (>= 3.10)
- R
- Bowtie2
- samtools
- MACS3 (Install below)

Ensure both Python and R are installed on your system.

## Installation

- R Dependencies
Install the required R packages by running the following commands in an R session:


```bash
git clone https://github.com/PM2951/ChIP-seq_Flow_Athaliana.git
cd ChIP-seq_AnalysisFlow
curl https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-60/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz -o Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
gunzip Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
bowtie2-build -f Arabidopsis_thaliana.TAIR10.dna.toplevel.fa /bowtie2_index/TAIR10
pip install openpyxl pandas customtkinter MACS3
```

### Rename directory name in config.py 

```Python
<your path> to your directory 
```
---

---

### parameter
- single end

  macs3 callpeak ; -p 0.05 -g 1.19e8

- pair end

  macs3 callpeak ; -p 0.05 -g 1.19e8 -f BAMPE

- annotation

  annotatePeakInBatch() ; featureType="TSS", multiple=TRUE, maxgap=10, PeakLocForDistance="middle", FeatureLocForDistance="TSS", select="first"

---

### License

This project is licensed under the Center for Gene Reserch, Nagoya Univ. (https://www.gene.nagoya-u.ac.jp/)

