# ChIP-seq_Flow_Athaliana

This repository provides a streamlined workflow for analyzing ChIP-seq data, leveraging Python and R tools to handle data preprocessing, visualization, and downstream analysis. Follow the instructions below to set up the environment and get started.

---

## Prerequisites

### Required Software
- Python (>= 3.8)
- R (>= 4.0)
- Bowtie2
- samtools
- MACS3

Ensure both Python and R are installed on your system.

## Installation

### Download files
```bash
curl https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-60/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz -o Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
gunzip Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
```
### Make index files

```bash
bowtie2-build -f Arabidopsis_thaliana.TAIR10.dna.toplevel.fa TAIR10
```

### Python Dependencies
Install the required Python packages using `pip`:

```bash
pip install openpyxl pandas
```

### R Dependencies
Install the required R packages by running the following commands in an R session:

```R
install.packages(c("dplyr", "ggplot2", "tidyr", "readxl"))
BiocManager::install("ChIPpeakAnno")
```

---

## How to Run the Workflow

1. Clone the repository:

   ```bash
   git clone https://github.com/PM2951/ChIP-seq_Flow_Athaliana.git
   echo "alias ChIPseq_main_run='python \"$HOME/PM2951/ChIPseq_main_run.py\"'" >> ~/.bashrc
   source ~/.bashrc
   ```

2. Execute the main script with required arguments.
  
   For example:
   
   ```bash
   ChIPseq_main_run -r 1 -n 32 -c Control_NAME -t Treated_NAME -m -p -a
   ```



   or (not alias ver.)

   ```bash
   python ChIPseq_main_run.py -r 1 -n 32 -c Control_NAME -t Treated_NAME -m -p -a
   ```


---

### License

This project is licensed under the MIT License. See the `LICENSE` file for more details.

