# ChIP-seq_Flow_Athaliana

This repository provides a streamlined workflow for analyzing ChIP-seq data, leveraging Python and R tools to handle data preprocessing, visualization, and downstream analysis. Follow the instructions below to set up the environment and get started.

---

## Prerequisites

### Required Software
- Python (>= 3.8)
- R (>= 4.0)
- Bowtie2
- samtools
- MACS3 (Install below)

Ensure both Python and R are installed on your system.

## Installation

- R Dependencies
Install the required R packages by running the following commands in an R session:

```R
install.packages(c("dplyr", "ggplot2", "tidyr", "readxl"))
BiocManager::install("ChIPpeakAnno")
```


```bash
git clone https://github.com/PM2951/ChIP-seq_Flow_Athaliana.git
cd ChIP-seq_AnalysisFlow
curl https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-60/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz -o Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
gunzip Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
bowtie2-build -f Arabidopsis_thaliana.TAIR10.dna.toplevel.fa TAIR10
pip install openpyxl pandas MACS3
```

---

## How to Run the Workflow

1. Clone the repository:

   ```bash
   git clone https://github.com/PM2951/ChIP-seq_Flow_Athaliana.git
   echo "alias ChIPseq_main_run='python \"$HOME/ChIP-seq_Flow_Athaliana/ChIPseq_main_run.py\"'" >> ~/.bashrc
   source ~/.bashrc
   ```

2. Move to your directory containing the required FASTQ files:
   ```bash
   cd <your directory>
   ```

   - FASTQ files are required (uncompressed .fastq format only).
    If your files are compressed (e.g., .gz), make sure to decompress them first using.
     ```bash
     gunzip *.gz
     ```

4. Execute the main script with required arguments.

   - help
     ```bash
     ChIPseq_main_run -h
     ```
  
   - main run:
   
      ```bash
      ChIPseq_main_run -r 1 -n 32 -c Control_NAME -t Treated_NAME -m -p -a
      ```



      or (not alias ver.)

      ```bash
      python ChIPseq_main_run.py -r 1 -n 32 -c Control_NAME -t Treated_NAME -m -p -a
      ```


   #### Command-line Options

   
   | Option               | Description                                                                                              |
   |----------------------|----------------------------------------------------------------------------------------------------------|
   | `-h`                 | Show the help message and exit.                                                                          |
   | `-r`                  | Specify read type: `1` for single-end reads, `2` for paired-end reads.                                   |
   | `-c`                  | Specify the common prefix (e.g., FILENAME) of the control sample files in FASTQ format. The files should follow the naming pattern: `FILENAME_001_R1_001.fastq`, `FILENAME_001_R1_002.fastq`, `FILENAME_001_R1_003.fastq`, `FILENAME_001_R1_004.fastq`. For example, if you input `FILENAME`, the script will process all matching files in the directory. The control sample is used to normalize the signal and reduce false positives in peak calling. |
   | `-t` 　　　　         | Specify the common prefix (e.g., FILENAME) of the treated sample files in FASTQ format. The files should follow the naming pattern: `FILENAME_001_R1_001.fastq`, `FILENAME_002_R1_001.fastq`, `FILENAME_003_R1_001.fastq`, `FILENAME_004_R1_001.fastq`. For example, if you input `FILENAME`, the script will process all matching files in the directory. The treated sample is used to identify regions with significant changes in signal compared to the control sample during peak calling. |
   | `-m`            　   | Perform read mapping to the reference genome.                                                           |
   | `-p`                 | Perform peak calling to identify ChIP-enriched regions.                                                 |
   | `-a`                 | Annotate identified peaks with genomic features (e.g., genes, promoters).                               |
   | `-n`                 | Specify number of CPU cores for parallel processing (default: 32).                                       |

   Replace the arguments as needed to suit your analysis.


---

### License

This project is licensed under the Center for Gene Reserch, Nagoya Univ. (https://www.gene.nagoya-u.ac.jp/)

