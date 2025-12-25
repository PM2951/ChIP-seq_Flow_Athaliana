import os

# config
# ChIP-seqアプリケーションの設定ファイル

# ディレクトリや実行ファイル
BOWTIE2_INDEX_DIR     = "Path your bowtie2 index directory"
BOWTIE2_BIN           = "bowtie2"
SAMTOOLS_BIN          = "samtools"
BEDTOOLS_BIN          = "bedtools"
MACS3_BIN             = "macs3"

# デフォルトのパラメータ
BOWTIE2_INDEX_NAME    = "TAIR10"
READ_TYPE             = "1"  # 1: single-end, 2: paired-end
DEFAULT_FOLDER_NAME   = "Path your favorite directory"

# リソースとパラメータ
NUM_CORES             = 15
DEFAULT_PVALUE        = "0.05"
GENOME_SIZE           = "1.19e8"

# MACS3 no modelの設定
MACS3_EXTSIZE         = 147

# ログ・スクリプト名
PARAMS_LOG_FILE       = "log_params.txt"
MAPPING_LOG_FILE      = "log_mapping.txt"
PEAKCALL_LOG_FILE     = "log_peakcall.txt"

# R アノテーションファイル
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
R_SCRIPT_NAME         = f"{SCRIPT_DIR}/R_annotation_code.R"
R_ANNOTATION_FILE     = f"{SCRIPT_DIR}/Athaliana_gene_ref.txt"