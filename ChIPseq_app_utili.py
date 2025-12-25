import subprocess
import pandas as pd
import os
import shutil
from pathlib import Path
from typing import List

import ChIPseq_app_config as config  

def check_executable(software_name: str):
    try:
        software_path = shutil.which(software_name)
        return software_path
    except Exception as e:
        raise RuntimeError(f"{software_name} not found. Install {software_name} or pass --{software_name} PATH")

def copy_if_not_exists(src_dir, dest_dir, filenames):
    for filename in filenames:
        src = os.path.join(src_dir, filename)
        dst = os.path.join(dest_dir, filename)
        if not os.path.exists(dst):
            shutil.copyfile(src, dst)


def write_run_parameters(output_dir, file_name=None, **kwargs):
    if file_name is None:
        file_name = config.PARAMS_LOG_FILE 
    param_log_file = os.path.join(output_dir, file_name)
    os.makedirs(output_dir, exist_ok=True)
    with open(param_log_file, "w", encoding="utf-8") as f:
        f.write("[ChIP-seq Run Parameters]\n")
        for key, value in kwargs.items():
            if isinstance(value, list):
                f.write(f"{key}:\n")
                for item in value:
                    f.write(f"  - {item}\n")
            else:
                f.write(f"{key}: {value}\n")


def mapping(output_name: str, input_files: List[str], output_dir: str, num_cores: int, index_name: str):
    # check software
    bowtie2 = check_executable(config.BOWTIE2_BIN)
    samtools = check_executable(config.SAMTOOLS_BIN)
    bedtools = check_executable(config.BEDTOOLS_BIN)

    r1_files = sorted([f for f in input_files if "_R1_" in Path(f).name])
    r2_files = sorted([f for f in input_files if "_R2_" in Path(f).name])

    if not r1_files:
        raise ValueError("R1ファイルが見つかりません")

    is_paired = bool(r2_files)

    r1_merged   = f"{output_dir}/{output_name}_R1.fastq"
    r2_merged   = f"{output_dir}/{output_name}_R2.fastq" if is_paired else None
    sam_file    = f"{output_dir}/{output_name}.sam"
    bam_file    = f"{output_dir}/{output_name}.bam"
    sorted_bam  = f"{output_dir}/{output_name}.sort.bam"


    if not os.path.exists(sam_file):
        with open(r1_merged, "wb") as w:
            for f in r1_files:
                with open(f, "rb") as r:
                    shutil.copyfileobj(r, w)

        if is_paired and r2_merged is not None:
            with open(r2_merged, "wb") as w:
                for f in r2_files:
                    with open(f, "rb") as r:
                        shutil.copyfileobj(r, w)

        MAPPING_INDEX = f"{config.BOWTIE2_INDEX_DIR}/{index_name}"
        cmd = [
            bowtie2,
            "-x", str(MAPPING_INDEX),
            "-S", str(sam_file),
            "-p", str(num_cores)
        ]
        
        if is_paired:
            cmd += ["-1", str(r1_merged), "-2", str(r2_merged)]
        else:
            cmd += ["-U", str(r1_merged)]

        outputlog_file = os.path.join(output_dir, config.MAPPING_LOG_FILE)
        
        mode = "a" if os.path.exists(outputlog_file) else "w"
        with open(outputlog_file, mode, encoding="utf-8") as log_file:
            result = subprocess.run(cmd, capture_output=True, text=True)
            log_file.write(f"Command: {' '.join(cmd)}\n")
            log_file.write(result.stdout)
            log_file.write(result.stderr)
            result.check_returncode()

            for samtools_cmd in [
                [samtools, "view", "-@", str(config.NUM_CORES), "-bS", str(sam_file), "-o", str(bam_file)],
                [samtools, "sort", "-@", str(config.NUM_CORES), str(bam_file), "-o", str(sorted_bam)],
                [samtools, "index", str(sorted_bam)],
                [bedtools, "genomecov", "-ibam", str(sorted_bam), "-bg", ">", f"{output_dir}/{output_name}.bedGraph"]
            ]:
                result = subprocess.run(samtools_cmd, capture_output=True, text=True)
                log_file.write(result.stdout)
                log_file.write(result.stderr)
                result.check_returncode()


def peakcall(output_dir, control_file_name, treated_file_name, read_type, pvalue, genome_size, no_model=False):
    # check software
    macs3 = check_executable(config.MACS3_BIN)

    if read_type not in ["1", "2"]:
        raise ValueError("read_type must be '1' (single-end) or '2' (paired-end)")

    peakcall_prefix = f"c{control_file_name}_t{treated_file_name}_peakcalling"
    peakcall_name   = f"{output_dir}/{peakcall_prefix}"
    peak_xls        = f"{output_dir}/{peakcall_prefix}_peaks.xls"
    peak_narrow     = f"{output_dir}/{peakcall_prefix}_peaks.narrowPeak"
    merged_xlsx     = f"{output_dir}/{peakcall_prefix}_peaks_merge.xlsx"

    if not os.path.exists(peak_xls):
        bam_control = f"{output_dir}/{control_file_name}.sort.bam"
        bam_treated = f"{output_dir}/{treated_file_name}.sort.bam"
        if not (os.path.exists(bam_control) and os.path.exists(bam_treated)):
            raise FileNotFoundError(f"BAM files not found: {bam_control}, {bam_treated}")

        cmd = [
            macs3, "callpeak",
            "-c", bam_control,
            "-t", bam_treated,
            "-n", peakcall_name,
            "-p", pvalue,
            "-g", genome_size,
        ]
        if read_type == "2":
            cmd += ["-f", "BAMPE"]

        if no_model:
            cmd += ["--nomodel", "--extsize", str(config.MACS3_EXTSIZE)]

        result = subprocess.run(cmd, capture_output=True, text=True)
        outputlog_file = os.path.join(output_dir, config.PEAKCALL_LOG_FILE)
        with open(outputlog_file, "a") as f:
            f.write(f"\n[MACS3 Output: {control_file_name} vs {treated_file_name}]\n")
            f.write(result.stdout)
            f.write(result.stderr)
        if result.returncode != 0:
            raise RuntimeError("MACS3 peak calling failed")

    if not os.path.exists(merged_xlsx):
        try:
            df_xls = pd.read_csv(peak_xls, sep='\t', skiprows=30)
            df_narrow = pd.read_table(peak_narrow, header=None)
            df_narrow.columns = [
                "chromosome", "start coordinate", "end coordinate", "name", "score",
                "strand", "signalValue", "pValue", "qValue", "peaksignal"
            ]
            df_merged = df_xls.merge(
                df_narrow[["name", "score", "peaksignal"]],
                on="name"
            )
            if "chr" in df_merged.columns:
                df_merged = df_merged[~df_merged["chr"].isin(["Mt", "Pt"])]
            df_merged.to_excel(merged_xlsx, index=False)
        except Exception as e:
            print(f"[ERROR] Excel merge failed: {e}")


def annotation(combined_excel_file):
    if not os.path.exists(f"{combined_excel_file}.xlsx"):
        raise FileNotFoundError(f"Combined Excel file not found: {combined_excel_file}")
    if not os.path.exists(config.R_SCRIPT_NAME):
        raise FileNotFoundError(f"R script not found: {config.R_SCRIPT_NAME}")
    if not os.path.exists(config.R_ANNOTATION_FILE):
        raise FileNotFoundError(f"R annotation file not found: {config.R_ANNOTATION_FILE}")
    
    subprocess.run(["chmod", "+x", config.R_SCRIPT_NAME], capture_output=True, text=True)

    command = ["Rscript", config.R_SCRIPT_NAME, combined_excel_file, config.R_ANNOTATION_FILE]
    try:
        result = subprocess.run(command, capture_output=True, check=True, text=True)
        print("Annotation succeeded.")
        print("\n", result.stdout)
    except subprocess.CalledProcessError as e:
        print("Error occurred during annotation:")
        print("Return code:", e.returncode)
        print("STDOUT:\n", e.stdout)
        print("STDERR:\n", e.stderr)


if __name__ == "__main__":
    print(
        """
control_file_name, treated_file_name: (example) YTHM21_S10, YTMN01_S3, ...
read_type: 1 or 2.  1 is single end, 2 is pair end.
mapping(control_file_name, treated_file_name, read_type) -> No return.
peakcall(control_file_name, treated_file_name, read_type) -> Returns Excel file name, used in annotation.
annotation(combined_excel_file) -> No return.
        """
    )
