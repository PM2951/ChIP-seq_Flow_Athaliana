print("ChIP-seq-appを起動中...")
print("10秒たっても何も表示されない場合は、一度消してから再度実行してください。")

# ChIP-seq解析のGUIアプリケーション
import os
import shutil
from time import sleep
import customtkinter as ctk
import tkinter.messagebox as messagebox
from tkinter import filedialog
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed

# 設定ファイルをインポート
import ChIPseq_app_config as config  
from ChIPseq_app_utili import mapping, peakcall, annotation, copy_if_not_exists


def write_run_parameters(output_dir, file_name=None, **kwargs):
    if file_name is None:
        file_name = config.PARAMS_LOG_FILE  # 遅延評価で安全
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


def askyesno_ctk(parent, title, message):
    dialog = ctk.CTkToplevel(parent)
    dialog.title(title)
    dialog.geometry("600x400")

    # --- ここを修正 ---
    dialog.update_idletasks()   # ウィンドウを表示準備
    dialog.grab_set()           # モーダル化
    # -------------------

    dialog.result = None

    frame = ctk.CTkFrame(dialog, corner_radius=12)
    frame.pack(expand=True, fill="both", padx=15, pady=15)

    title_label = ctk.CTkLabel(
        frame, text=title,
        font=ctk.CTkFont(size=18, weight="bold"),
        text_color="#e63946"
    )
    title_label.pack(pady=(10, 5))

    msg_label = ctk.CTkLabel(
        frame, text=message,
        wraplength=550, justify="center",
        font=ctk.CTkFont(size=14)
    )
    msg_label.pack(pady=10, padx=10)

    btn_frame = ctk.CTkFrame(frame, fg_color="transparent")
    btn_frame.pack(pady=10)

    def yes():
        dialog.result = True
        dialog.destroy()

    def no():
        dialog.result = False
        dialog.destroy()

    ctk.CTkButton(btn_frame, text="実行する", width=100,
                  fg_color="#2a9d8f", hover_color="#21867a",
                  command=yes).pack(side="left", padx=10)

    ctk.CTkButton(btn_frame, text="キャンセル", width=100,
                  fg_color="#e76f51", hover_color="#d65a3a",
                  command=no).pack(side="left", padx=10)

    parent.wait_window(dialog)
    return dialog.result


# GUI部分
class App(ctk.CTk):
    def __init__(self):
        super().__init__()
        self.title("ChIP-seq All-in-One")
        self.geometry("1200x800")

        # 状態変数
        self.control_files = []
        self.treatment_files = []
        self.output_dir = ""  # 出力フォルダパス

        # GUI変数
        self.control_name_var   = ctk.StringVar()
        self.treatment_name_var = ctk.StringVar()  # デフォルト値を設定
        self.mapping_var        = ctk.BooleanVar(value=True)
        self.peakcall_var       = ctk.BooleanVar(value=True)
        self.annotation_var     = ctk.BooleanVar(value=True)
        self.no_model_var      = ctk.BooleanVar(value=False)
        
        self.read_type_var      = ctk.StringVar(value=config.READ_TYPE)
        self.index_name_var     = ctk.StringVar(value=config.BOWTIE2_INDEX_NAME)
        self.pvalue_var         = ctk.StringVar(value=config.DEFAULT_PVALUE)  
        self.genome_size_var    = ctk.StringVar(value=config.GENOME_SIZE) 

        self._create_widgets()


    def _create_widgets(self):
        # 上部フレーム（左右に分ける）
        top_frame = ctk.CTkFrame(self)
        top_frame.pack(side="top", fill="both", expand=True, padx=10, pady=10)

        top_frame.grid_columnconfigure(0, weight=2)  
        top_frame.grid_columnconfigure(1, weight=8) 
        top_frame.grid_rowconfigure(0, weight=1)

        left_frame = ctk.CTkFrame(top_frame)
        left_frame.grid(row=0, column=0, sticky="nsew", padx=(0, 5))  # 左右パディング

        right_frame = ctk.CTkFrame(top_frame)
        right_frame.grid(row=0, column=1, sticky="nsew", padx=(5, 0))

        # right
        Control_name_frame = ctk.CTkFrame(right_frame)
        Control_name_frame.pack(pady=(10, 10))
        ctk.CTkLabel(Control_name_frame, text="Controlサンプル名").pack(side="left", padx=(10, 10))
        ctk.CTkEntry(Control_name_frame, textvariable=self.control_name_var).pack(side="left")
        ctk.CTkButton(right_frame, text="Select Control Files", command=self.select_control).pack(pady=(10, 10))
        self.control_listbox = ctk.CTkTextbox(right_frame, height=60)
        self.control_listbox.insert("end", "PeakcallまたはAnnotationの場合は、sort.bamのファイル名を入力してください。\n参考) FILE_NAME.sort.bam -> FILE_NAME\n")
        self.control_listbox.pack(pady=(10,10), padx=(25,25), fill="both", expand=True)

        treatment_name_frame = ctk.CTkFrame(right_frame)
        treatment_name_frame.pack(pady=(10, 10))
        ctk.CTkLabel(treatment_name_frame, text="treatmentサンプル名").pack(side="left", padx=(10, 10))
        ctk.CTkEntry(treatment_name_frame, textvariable=self.treatment_name_var).pack(side="left")
        ctk.CTkButton(right_frame, text="Select Treatment Files", command=self.select_treatment).pack(pady=(10, 10))
        self.treatment_listbox = ctk.CTkTextbox(right_frame, height=60)
        self.treatment_listbox.insert("end", "PeakcallまたはAnnotationの場合は、sort.bamのファイル名を入力してください。\n参考) FILE_NAME.sort.bam -> FILE_NAME\n")
        self.treatment_listbox.pack(pady=(10,10), padx=(25,25), fill="both", expand=True)

        # オプションや出力関連
        option_frame = ctk.CTkFrame(right_frame)
        option_frame.pack(fill="x", padx=25, pady=10)
        ctk.CTkButton(option_frame, text="出力フォルダを選択", command=self.select_output_dir).pack(pady=(30, 0))
        self.output_dir_label = ctk.CTkLabel(option_frame, text="未選択", wraplength=700)
        self.output_dir_label.pack(padx=10, pady=15)


        # left
        # 中央フレーム（オプション関連）
        ctk.CTkCheckBox(left_frame, text="Mapping", variable=self.mapping_var).pack(anchor="center", pady=(40, 0))
        ctk.CTkCheckBox(left_frame, text="Peak Calling", variable=self.peakcall_var).pack(anchor="center", pady=2)
        ctk.CTkCheckBox(left_frame, text="Annotation", variable=self.annotation_var).pack(anchor="center", pady=2)

        # Mapping options
        ctk.CTkLabel(left_frame, text="Mapping options").pack(pady=(20, 0), anchor="center")

        radio_frame = ctk.CTkFrame(left_frame)
        radio_frame.pack(pady=5)
        ctk.CTkRadioButton(radio_frame, text="Single-end", variable=self.read_type_var, value="1").pack(side="left", padx=10)
        ctk.CTkRadioButton(radio_frame, text="Paired-end", variable=self.read_type_var, value="2").pack(side="left", padx=10)
        
        # Mapping Index Name 入力欄
        index_frame = ctk.CTkFrame(left_frame)
        index_frame.pack(pady=(10, 0))
        ctk.CTkLabel(index_frame, text="Index（Bowtie2）:").pack(side="left", padx=(5, 5))
        ctk.CTkEntry(index_frame, textvariable=self.index_name_var, width=120).pack(side="left")
        # ctk.CTkLabel(left_frame, text="indexファイルは、\n「pm2951/ChIP-seq_Flow_Athaliana/bowtie2_index」\nに入れてください。").pack(pady=(0, 0), anchor="center")

        # peakcall p-value
        ctk.CTkLabel(left_frame, text="Peakcall options").pack(pady=(20, 0), anchor="center")

        pval_frame = ctk.CTkFrame(left_frame)
        pval_frame.pack(pady=(5, 0))
        ctk.CTkLabel(pval_frame, text="p-value:").pack(side="left", padx=(5, 5))
        ctk.CTkEntry(pval_frame, textvariable=self.pvalue_var, width=80).pack(side="left", padx=10)

        # peakcall genome size
        genome_frame = ctk.CTkFrame(left_frame)
        genome_frame.pack(pady=(5, 10))
        ctk.CTkLabel(genome_frame, text="Genomeサイズ (1.23e4など):").pack(side="left", padx=(5, 5))
        ctk.CTkEntry(genome_frame, textvariable=self.genome_size_var, width=80).pack(side="left", padx=10)

        # no model
        ctk.CTkCheckBox(left_frame, text="MACS3 no model (if needed)", variable=self.no_model_var).pack(anchor="center", pady=2)

        # 実行ボタン
        ctk.CTkButton(self, text="解析を実行", command=self.run_process).pack(pady=10)

        # 実行ログ（ウィンドウ下部）
        self.progressbar = ctk.CTkProgressBar(self, 
                                            height=10,      # 高さ
                                            width=300,      # 幅
                                            )

        self.progressbar.pack_forget()  # 最初は非表示にする
    

    def select_control(self):
        start_dir = os.path.expanduser(config.DEFAULT_FOLDER_NAME)
        
        files = filedialog.askopenfilenames(
            title="Control fastaファイル選択", 
            initialdir=start_dir,  # 初期フォルダ指定
            filetypes=[("FASTA/FASTQ Files", "*.fa *.fasta *.fq *.fastq *.gz")])
        if files:
            self.control_files = list(files)
            self.control_listbox.delete("0.0", "end")
            for f in self.control_files:
                self.control_listbox.insert("end", f"{os.path.basename(f)}\n")

    def select_treatment(self):
        start_dir = os.path.expanduser(config.DEFAULT_FOLDER_NAME)

        files = filedialog.askopenfilenames(
            title="Treatment fastaファイル選択", 
            initialdir=start_dir,  # 初期フォルダ指定
            filetypes=[("FASTA/FASTQ Files", "*.fa *.fasta *.fq *.fastq *.gz")])
        if files:
            self.treatment_files = list(files)
            self.treatment_listbox.delete("0.0", "end")
            for f in self.treatment_files:
                self.treatment_listbox.insert("end", f"{os.path.basename(f)}\n")


    def select_output_dir(self):
        start_dir = os.path.expanduser(config.DEFAULT_FOLDER_NAME)

        directory = filedialog.askdirectory(
            title="出力フォルダを選択",
            initialdir=start_dir  # 初期フォルダ指定
        )
        if directory:
            self.output_dir = directory
            self.output_dir_label.configure(text=self.output_dir)

    def run_process(self):

        # サンプル名が未入力 → エラーダイアログ（showerror）
        if not self.control_name_var.get() or not self.treatment_name_var.get():
            messagebox.showerror(
                "サンプル名エラー",
                "ControlまたはTreatmentのサンプル名が未入力です。",
                parent=self
            )
            return
        
        # Control/Treatmentファイルが未選択 → 確認ダイアログ（askyesno）
        if not self.control_files or not self.treatment_files:
            proceed = askyesno_ctk(
                self,
                "ファイル未選択\n",
                "ControlまたはTreatmentのファイルが選択されていません。"
                "PeakCallまたはAnnotationは未選択でも実行できます。\n\n"
                "！！！注意！！！\n\n"
                "PeakcallまたはAnnotationの場合は、sort.bamのファイル名を入力してください。\n"
                "参考) FILE_NAME.sort.bam -> FILE_NAME\n\n\n"
                "このまま実行しますか？"
            )

            if proceed and self.mapping_var.get():
                messagebox.showerror(
                    "mappingにチェックが入っています",
                    "チェックを外すか、ControlまたはTreatmentのファイルを選択してください。",
                    parent=self
                )
                return

            elif not proceed:
                return

        
        
        self.progressbar.pack(pady=40)
        self.progressbar.set(0)

        from datetime import datetime, timedelta
        finish_time = datetime.now() + timedelta(minutes=10)
        print("\n完了予定時刻:", finish_time.strftime("%Y-%m-%d %H:%M"), flush=True)

        print("解析ファイル：", flush=True)
        print(" -Control files:")
        for i, f in enumerate(self.control_files):
            print(f"    - {i+1} | {f}")

        print(" -Treatment files:")
        for i, f in enumerate(self.treatment_files):
            print(f"    - {i+1} | {f}")
        print(" -Output directory:")
        print(f"    - {self.output_dir}", flush=True)

        print("\n解析を実行中.....", flush=True)

        self.progressbar.set(0.2)


        # 解析を実行
        params = {
            "mapping_flag": self.mapping_var.get(),
            "peakcall_flag": self.peakcall_var.get(),
            "annotation_flag": self.annotation_var.get(),
            "control_name": self.control_name_var.get(),
            "treated_name": self.treatment_name_var.get(),
            "control_files": self.control_files,
            "treatment_files": self.treatment_files,
            "output_dir": self.output_dir,
            "read_type": str(self.read_type_var.get()),
            "index_name": self.index_name_var.get(),
            "pvalue": self.pvalue_var.get(),
            "genome_size": self.genome_size_var.get(),
            "num_cores": config.NUM_CORES,
            "no_model": self.no_model_var.get(),
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

            self.progressbar.set(0.6)

            if params["peakcall_flag"]:
                peakcall(params["output_dir"], params["control_name"], params["treated_name"], params["read_type"],
                        pvalue=params["pvalue"], genome_size=params["genome_size"],
                        no_model=params["no_model"])
                print("[Peak Calling] Peak calling finished.", flush=True)
            self.progressbar.set(0.8)

            if params["annotation_flag"]:
                peak_file = os.path.join(params["output_dir"], f"c{params['control_name']}_t{params['treated_name']}_peakcalling_peaks_merge")    #Rの中で.xlsxに変換される
                annotation(peak_file)
                print("[Annotation] Annotation finished.", flush=True)
            self.progressbar.set(0.9)
            print("\nAll processes completed successfully!", flush=True)
        except Exception as e:
            messagebox.showerror(
                "解析エラー",
                f"解析中にエラーが発生しました:\n{e}",
                parent=self
            )
            print(f"Error during processing: {e}", flush=True)
            self.progressbar.set(0.0)
            return


        print("\n=== Finished ===", flush=True)
        self.progressbar.set(1.0)

        sleep(1)
        self.progressbar.pack_forget()  # 最初は非表示にする

if __name__ == "__main__":
    ctk.set_appearance_mode("light")
    app = App()
    print("=== Welcome to ChIP-seq-app ===")


    app.mainloop()
