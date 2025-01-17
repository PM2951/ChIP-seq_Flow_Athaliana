
import subprocess
import pandas as pd
import os
import shutil

current_directory = os.getcwd()
files_in_directory = os.listdir(current_directory)

script_path = os.path.dirname(os.path.abspath(__file__))
singe_script_path = f"{script_path}/mapping_single.sh"
pairend_script_path = f"{script_path}/mapping_pairend.sh"
r_script_path = f"{script_path}/R_annotation_code.R"

# Grant execution permissions to scripts
subprocess.run(["chmod", "+x", singe_script_path], capture_output=True, text=True)
subprocess.run(["chmod", "+x", pairend_script_path], capture_output=True, text=True)
subprocess.run(["chmod", "+x", r_script_path], capture_output=True, text=True)

def print_green(text):
    """Prints text in green color.

    Args:
        text (str): The text to print.
    """
    print("\033[92m" + text + "\033[0m")

def print_red(text):
    """Prints text in red color.

    Args:
        text (str): The text to print.
    """
    print("\033[91m" + text + "\033[0m")


def mapping(fastaq_file_name, read_type, num_cores):
    """Performs read mapping using a shell script.

    Args:
        fastaq_file_name (str): Name of the input FASTQ file (without extension).
        read_type (str): '1' for single-end reads, '2' for paired-end reads.

    Raises:
        TypeError: If read_type is invalid or mapping fails.
    """
    print_green("========== mapping ==========")
    mapping_result = None
    try:
        if f"{fastaq_file_name}.sam" not in files_in_directory:
            if read_type == "1":
                script = singe_script_path
            elif read_type == "2":
                script = pairend_script_path
            else:
                raise TypeError('Invalid read_type argument.')

            mapping_result = subprocess.run([script, fastaq_file_name, num_cores],
                                             capture_output=True, text=True)

            output_file = 'output_mapping.txt'
            subprocess.run(["chmod", "+w", output_file], capture_output=True, text=True)
            mode = "a" if os.path.exists(output_file) else "w"
            with open(output_file, mode) as f:
                f.write(f"\n{fastaq_file_name}\n")
                f.write(mapping_result.stdout)
                f.write(mapping_result.stderr)
                f.write("\n")
            print("Mapping results saved in", output_file)
    except:
        raise TypeError('Mapping error')

def peakcall(contol_file_name, treated_file_name, read_type):
    """Performs peak calling using MACS3.

    Args:
        contol_file_name (str): Control file base name (without extension).
        treated_file_name (str): Treated file base name (without extension).
        read_type (str): '1' for single-end reads, '2' for paired-end reads.

    Raises:
        TypeError: If read_type is invalid or peak calling fails.
    """
    print_green("========== peakcalling ==========")
    if f"c{contol_file_name}_t{treated_file_name}_peakcalling_peaks.xls" not in files_in_directory:
        peakcall_result = None
        if read_type == "1":
            peakcall_result = subprocess.run(["macs3", "callpeak",
                                              "-c", f"{contol_file_name}.sort.bam",
                                              "-t", f"{treated_file_name}.sort.bam",
                                              "-n", f"c{contol_file_name}_t{treated_file_name}_peakcalling",
                                              "-p", "0.05",
                                              "-g", "1.19e8"],
                                             capture_output=True, text=True)
        elif read_type == "2":
            peakcall_result = subprocess.run(["macs3", "callpeak",
                                              "-c", f"{contol_file_name}.sort.bam",
                                              "-t", f"{treated_file_name}.sort.bam",
                                              "-n", f"c{contol_file_name}_t{treated_file_name}_peakcalling",
                                              "-p", "0.05",
                                              "-g", "1.19e8",
                                              "-f", "BAMPE"],
                                             capture_output=True, text=True)
        else:
            raise TypeError('Peakcalling error')

        output_file = 'output_peakcall.txt'
        subprocess.run(["chmod", "+w", output_file], capture_output=True, text=True)
        mode = "a" if os.path.exists(output_file) else "w"
        with open(output_file, mode) as f:
            f.write(f"\ncontrol: {contol_file_name} \ntreatment: {treated_file_name}\n")
            f.write(peakcall_result.stdout)
            f.write(peakcall_result.stderr)
        print("Peak calling results saved in", output_file)

    if f"c{contol_file_name}_t{treated_file_name}_peakcalling_peaks_merge.xlsx" not in files_in_directory:
        try:
            print(" ==> Merging Excel files...")
            file = f"c{contol_file_name}_t{treated_file_name}_peakcalling_peaks.xls"
            df = pd.read_csv(file, sep='\t', skiprows=30)
            file_narrow = f"c{contol_file_name}_t{treated_file_name}_peakcalling_peaks.narrowPeak"
            df_narrow = pd.read_table(file_narrow, header=None)
            df_narrow.columns = ["chromosome", "start coordinate", "end coordinate", "name", "score", "strand", "signalValue", "pValue", "qValue", "peaksignal"]
            df_merge = df.merge(df_narrow[['score', 'peaksignal']], left_index=True, right_index=True)
            df_merge = df_merge[(df_merge['chr'] != 'Mt') & (df_merge['chr'] != 'Pt')]
            outfile = f"c{contol_file_name}_t{treated_file_name}_peakcalling_peaks_merge.xlsx"
            df_merge.to_excel(outfile, index=False)
            print(f"Saved merged file as {outfile}")
        except:
            print_red("Error occurred during Excel processing")

def annotation(combined_excel_file):
    """Performs annotation using an R script.

    Args:
        combined_excel_file (str): Path to the merged Excel file from peak calling.

    Raises:
        subprocess.CalledProcessError: If the R script execution fails.
    """
    print_green("========== annotation ==========")
    arguments = [combined_excel_file]
    command = ["Rscript", r_script_path] + arguments
    try:
        result = subprocess.run(command, capture_output=True, check=True, text=True)
        print("Annotation completed successfully")
        print(result.stdout)
    except subprocess.CalledProcessError as e:
        print_red("Error occurred during annotation:")
        print(e.stderr)

if __name__ == "__main__":
    print_green(
        """
control_file_name, treated_file_name: (example) YTHM21_S10, YTMN01_S3, ...
read_type: 1 or 2.  1 is single end, 2 is pair end.
mapping(control_file_name, treated_file_name, read_type) -> No return.
peakcall(control_file_name, treated_file_name, read_type) -> Returns Excel file name, used in annotation.
annotation(combined_excel_file) -> No return.
        """
    )
