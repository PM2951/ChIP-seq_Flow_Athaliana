import os
import shutil
from time import sleep
from datetime import datetime, timedelta
from concurrent.futures import ThreadPoolExecutor, as_completed
import argparse

# 設定ファイルをインポート
import config  
from ChIPseq_app_utili import mapping, peakcall, annotation, write_run_parameters


file_extentions = (".fastq", ".fastq.gz", ".fq", ".fq.gz", "fasta", ".fasta.gz")


# Create parser with a detailed description of the pipeline
p = argparse.ArgumentParser(
    description=(
        "ChIP-seq analysis pipeline to process sequencing data, "
        "including mapping, peak calling, and annotation. "
        "Specify input files and choose the analysis steps to perform."
    ))

p.add_argument(
    "-r", '--read_type', choices=['1', '2'], default=config.READ_TYPE,
    help=(
        "Specify the type of sequencing reads to process. "
        "Choose '1' for single-end reads or '2' for paired-end reads. "
        "This affects the alignment and downstream analysis steps."
    ))

p.add_argument(
    "-c", '--control_name', required=True,
    help=(
        "Specify the common prefix (e.g., FILENAME) of the control sample files in FASTQ format. "
        "The files should follow the naming pattern: "
        "'FILENAME_001_R1_001.fastq', 'FILENAME_002_R1_001.fastq', 'FILENAME_001_R1_003.fastq', 'FILENAME_002_R1_004.fastq'. "
        "For example, if you input 'FILENAME', the script will process all matching files in the directory. "
        "The control sample is used to normalize the signal and reduce false positives in peak calling."
        f"supported file extensions are: {', '.join(file_extentions)}"
    ))

p.add_argument(
    "-t", '--treated_name', required=True,
    help=(
        "Specify the common prefix (e.g., FILENAME) of the treated sample files in FASTQ format."
        "The files should follow the naming pattern:"
        "  'FILENAME_001_R1_001.fastq', 'FILENAME_002_R1_001.fastq', 'FILENAME_003_R1_001.fastq', 'FILENAME_004_R1_001.fastq'."
        "For example, if you input 'FILENAME', the script will process all matching files in the directory."
        "The treated sample is used to identify regions with significant changes in signal compared to the control sample during peak calling."
        f"supported file extensions are: {', '.join(file_extentions)}"
    ))

p.add_argument(
    "-o", '--output_dir', default="output",
    help=(
        "Specify the output directory for the analysis results. "
        "All output files will be saved in this directory."
    ))

p.add_argument(
    "-m", '--mapping_flag', action='store_true',
    help=(
        "Flag to perform read mapping to the reference genome. "
        "This step aligns sequencing reads to the reference and outputs BAM files."
    ))

p.add_argument(
    "-p", '--peakcall_flag', action='store_true',
    help=(
        "Flag to perform peak calling. "
        "This step identifies regions of significant ChIP enrichment in the genome."
    ))

p.add_argument(
    "-a", '--annotation_flag', action='store_true',
    help=(
        "Flag to perform annotation of the identified peaks. "
        "The peaks will be annotated with genomic features (e.g., genes, promoters)."
    ))
p.add_argument(
    "-n", '--num_cores', type=int, default=config.NUM_CORES,
    help=(
        "Number of CPU cores to use for parallel processing. "
        "Specify this to speed up the analysis by using multiple cores."
        f"Default is {config.NUM_CORES}."
    ))
p.add_argument(
    "-i", '--index_name', default=config.BOWTIE2_INDEX_NAME,
    help=(
        "Specify the name of the Bowtie2 index to use for read mapping. "
        "This should correspond to the reference genome you are aligning against. "
        f"Default is {config.BOWTIE2_INDEX_NAME} for Arabidopsis thaliana."
    ))
p.add_argument(
    "-g", '--genome_size', default=config.GENOME_SIZE,
    help=(
        "Effective genome size for peak calling. "
        "This is used by MACS3 to estimate the background model. "
        f"Default is {config.GENOME_SIZE} for Arabidopsis thaliana."
    ))
p.add_argument(
    "-pv", '--pvalue', default=config.DEFAULT_PVALUE,
    help=(
        "P-value threshold for peak calling. "
        "This sets the stringency for identifying significant peaks. "
        f"Default is {config.DEFAULT_PVALUE}"
    ))
p.add_argument(
    "--no_model", action='store_true',
    help=(
        "Flag to skip model building during peak calling. "
        "Use this if you want to call peaks without the default model."
    ))

args = p.parse_args()

current_dir = os.getcwd()
output_dir = os.path.join(current_dir, args.output_dir)
os.makedirs(output_dir, exist_ok=True)


# コマンドライン引数の表示
print(" Parameters")
print(f"    -read_type  : {args.read_type}")
print(f"    -control    : {args.control_name}")
print(f"    -treated    : {args.treated_name}")
print(f"    -output     : {args.output_dir}")
print(f"    -mapping    : {args.mapping_flag}")
print(f"    -peakcall   : {args.peakcall_flag}")
print(f"    -annotation : {args.annotation_flag}")
print(f"    -core       : {args.num_cores}")


control_files = sorted([f for f in os.listdir(current_dir) if f.startswith(args.control_name) and f.endswith(file_extentions)])
if len(control_files) == 0:
    raise FileNotFoundError(f"Error: No control FASTQ files found with prefix '{args.control_name}'")

treatment_files = sorted([f for f in os.listdir(current_dir) if f.startswith(args.treated_name) and f.endswith(file_extentions)])
if len(treatment_files) == 0:
    raise FileNotFoundError(f"Error: No treatment FASTQ files found with prefix '{args.treated_name}'")

print("FILES:", flush=True)
print(" -Control files:")
for i, f in enumerate(control_files):
    print(f"    - {i+1}: {f}")

print(" -Treatment files:")
for i, f in enumerate(treatment_files):
    print(f"    - {i+1}: {f}")

finish_time = datetime.now() + timedelta(minutes=10)
print("\nExpected finish time:", finish_time.strftime("%Y-%m-%d %H:%M"), flush=True)
print("\nRunning analysis.....", flush=True)



# 解析を実行
params = {
    "mapping_flag": args.mapping_flag,
    "peakcall_flag": args.peakcall_flag,
    "annotation_flag": args.annotation_flag,
    "control_name": args.control_name,
    "treated_name": args.treated_name,
    "control_files": control_files,
    "treatment_files": treatment_files,
    "output_dir": output_dir,
    "read_type": str(args.read_type),
    "index_name": args.index_name,
    "pvalue": args.pvalue,
    "genome_size": args.genome_size,
    "num_cores": config.NUM_CORES,
    "no_model": args.no_model,
}

# output_dirを除いたパラメータ辞書
params_without_output_dir = {k: v for k, v in params.items() if k != "output_dir"}
write_run_parameters(
    output_dir=params["output_dir"],
    file_name=f"c{params['control_name']}_t{params['treated_name']}_{config.PARAMS_LOG_FILE}",
    **params_without_output_dir
)

try:
    if params["mapping_flag"]:
        print("[Mapping] Running control and treatment mapping in parallel...", flush=True)

        with ThreadPoolExecutor(max_workers=2) as executor:
            futures = [
                executor.submit(mapping, params["control_name"], params["control_files"], params["output_dir"], params["num_cores"], params["index_name"]),
                executor.submit(mapping, params["treated_name"], params["treatment_files"], params["output_dir"], params["num_cores"], params["index_name"]),
            ]

            # Ensure that any exceptions in the worker threads are raised here
            for f in as_completed(futures):
                f.result()
        print("[Mapping] Both mapping jobs finished.", flush=True)
except Exception as e:
    raise ValueError(f"Error during processing: {e}")

try:
    if params["peakcall_flag"]:
        peakcall(params["output_dir"], params["control_name"], params["treated_name"], params["read_type"],
                pvalue=params["pvalue"], genome_size=params["genome_size"],
                no_model=params["no_model"])
        print("[Peak Calling] Peak calling finished.", flush=True)
except:
    log_file_path = os.path.join(params["output_dir"], "log_peakcall.txt")
    if os.path.exists(log_file_path):
        with open(log_file_path, "r") as log_file:
            lines = log_file.readlines()
            print("".join(lines[-3:]), flush=True)
    raise ValueError("Error during peak calling. Check log_peakcall.txt for details.")

try:
    if params["annotation_flag"]:
        peak_file = os.path.join(params["output_dir"], f"c{params['control_name']}_t{params['treated_name']}_peakcalling_peaks_merge")    #Rの中で.xlsxに変換される
        annotation(peak_file)
        print("[Annotation] Annotation finished.", flush=True)

    print("\nAll processes completed successfully!", flush=True)
except Exception as e:
    raise ValueError(f"Error during processing: {e}")


sleep(1)
