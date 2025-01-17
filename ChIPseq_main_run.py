#解析したfastaqがあるファルダにcdで移動
#MN05_S18などの4分割されたfastaqに共通して存在するファイル名を入力してrun

import subprocess
import os
import shutil
from ChIPseq_utili import mapping, peakcall, annotation
import argparse

#コマンドラインに色をつける
def print_green(text):
    print("\033[92m" + text + "\033[0m")
def print_red(text):
    print("\033[91m" + text + "\033[0m")
def copy_if_not_exists(src_dir, dest_dir, filenames):
    for filename in filenames:
        src_file = os.path.join(src_dir, filename)
        dest_file = os.path.join(dest_dir, filename)
        if not os.path.exists(dest_file):
            shutil.copyfile(src_file, dest_file)
            print(f"File '{filename}' copied")
def making_dir(output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"Output directory '{output_dir}' created")


# Create parser with a detailed description of the pipeline
parser = argparse.ArgumentParser(
    description=(
        "ChIP-seq analysis pipeline to process sequencing data, "
        "including mapping, peak calling, and annotation. "
        "Specify input files and choose the analysis steps to perform."
    ))
# Add arguments with detailed help messages
parser.add_argument(
    "-r", '--read_type', choices=['1', '2'], required=True,
    help=(
        "Specify the type of sequencing reads to process. "
        "Choose '1' for single-end reads or '2' for paired-end reads. "
        "This affects the alignment and downstream analysis steps."
    ))
parser.add_argument(
    "-c", '--control', required=True,
    help=(
    "Specify the common prefix (e.g., FILENAME) of the control sample files in FASTQ format. "
    "The files should follow the naming pattern: "
    "'FILENAME_001_R1_001.fastq', 'FILENAME_002_R1_001.fastq', 'FILENAME_001_R1_003.fastq', 'FILENAME_002_R1_004.fastq'. "
    "For example, if you input 'FILENAME', the script will process all matching files in the directory. "
    "The control sample is used to normalize the signal and reduce false positives in peak calling."
    ))
parser.add_argument(
    "-t", '--treated', required=True,
    help=(
        "Specify the common prefix (e.g., FILENAME) of the treated sample files in FASTQ format."
        "The files should follow the naming pattern:"
        "  'FILENAME_001_R1_001.fastq', 'FILENAME_002_R1_001.fastq', 'FILENAME_003_R1_001.fastq', 'FILENAME_004_R1_001.fastq'."
        "For example, if you input 'FILENAME', the script will process all matching files in the directory."
        "The treated sample is used to identify regions with significant changes in signal compared to the control sample during peak calling."
    ))
parser.add_argument(
    "-m", '--mapping', action='store_true',
    help=(
        "Flag to perform read mapping to the reference genome. "
        "This step aligns sequencing reads to the reference and outputs BAM files."
    ))

parser.add_argument(
    "-p", '--peakcall', action='store_true',
    help=(
        "Flag to perform peak calling. "
        "This step identifies regions of significant ChIP enrichment in the genome."
    ))

parser.add_argument(
    "-a", '--annotation', action='store_true',
    help=(
        "Flag to perform annotation of the identified peaks. "
        "The peaks will be annotated with genomic features (e.g., genes, promoters)."
    ))
parser.add_argument(
    "-n", '--num_cores', type=int, default=32,
    help=(
        "Number of CPU cores to use for parallel processing. "
        "Specify this to speed up the analysis by using multiple cores."
        "default 32"
    ))
args = parser.parse_args()

# コマンドライン引数の表示
print(" Parameters")
print(f"    -read_type  : {args.read_type}")
print(f"    -control    : {args.control}")
print(f"    -treated    : {args.treated}")
print(f"    -mapping    : {args.mapping}")
print(f"    -peakcall   : {args.peakcall}")
print(f"    -annotation : {args.annotation}")
print(f"    -core       : {args.num_cores}")

# 出力ディレクトリの作成
# ディレクトリとファイル名の設定, ChIP用TAIR10ファイルの取得
src_dir = os.path.dirname(os.path.abspath(__file__))  # ChIP用TAIR10ファイルのディレクトリ
current_dir = os.getcwd()
filenames = ["Arabidopsis_thaliana.TAIR10.dna.toplevel.fa", 
             "Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.fai", 
             "TAIR10.1.bt2", "TAIR10.2.bt2", "TAIR10.3.bt2", "TAIR10.4.bt2", 
             "TAIR10.rev.1.bt2", "TAIR10.rev.2.bt2",
             "Athaliana_gene_ref.txt",
             ]
copy_if_not_exists(src_dir, current_dir, filenames)

# # main functions
if args.mapping:
    mapping(args.control, args.read_type, str(args.num_cores))
    mapping(args.treated, args.read_type, str(args.num_cores))

if args.peakcall:
    peakcall(args.control, args.treated, args.read_type)

if args.annotation:
    annotation(f"{current_dir}/c{args.control}_t{args.treated}_peakcalling_peaks_merge")






# # one time run, ファイルが大きすぎると途中でbowtie2が処理落ちする可能性あり。
# read_type = "1"
# file_names = [
#     "YTMN14_S2",
# ]
# for file_name in file_names:
#     control_file_name = "YTMN16_S4"
#     treated_file_name = file_name
#     peakcall(control_file_name, treated_file_name, read_type)
#     annotation(f"{current_dir}/c{control_file_name}_t{treated_file_name}_peakcalling_peaks_merge")