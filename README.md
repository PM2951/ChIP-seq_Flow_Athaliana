# ChIP-seq_Flow_Athaliana

This repository provides a streamlined workflow for analyzing ChIP-seq data, leveraging Python and R tools to handle data preprocessing, visualization, and downstream analysis. Follow the instructions below to set up the environment and get started.

---

## Prerequisites

### Required Software
- Python (>= 3.8)
- R (>= 4.0)

Ensure both Python and R are installed on your system.

## Installation

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
   git clone https://github.com/<username>/ChIP-seq_AnalysisFlow.git
   echo "alias ChIPseq_main_run='python \"$HOME/repository-name/ChIPseq_main_run.py\"'" >> ~/.bashrc
   source ~/.bashrc
   ```

2. Execute the main script with required arguments.
  
   For example:

   ```bash
   python ChIPseq_main_run.py -r 1 -n 32 -c Control_NAME -t Treated_NAME -m -p -a
   ```


---

### License

This project is licensed under the MIT License. See the `LICENSE` file for more details.

