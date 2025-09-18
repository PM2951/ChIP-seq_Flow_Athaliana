# config
# ChIP-seqアプリケーションの設定ファイル

# ディレクトリや実行ファイル
BOWTIE2_INDEX_DIR     = "/home/pm2951/ChIP-seq_Flow_Athaliana/bowtie2_index"
BOWTIE2_BIN           = "bowtie2"
SAMTOOLS_BIN          = "samtools"
MACS3_BIN             = "macs3"

# デフォルトのパラメータ
BOWTIE2_INDEX_NAME    = "TAIR10"
READ_TYPE             = "1"  # "single-end" または "paired-end"     
DEFAULT_FOLDER_NAME   = "/home/pm2951/shortcut"

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
R_SCRIPT_NAME         = "/home/pm2951/ChIP-seq_Flow_Athaliana/R_annotation_code.R"
R_ANNOTATION_FILE     = "/home/pm2951/ChIP-seq_Flow_Athaliana/Athaliana_gene_ref.txt"
