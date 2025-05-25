from handler.utils import *
from handler.reader import *
from handler.config_parser import get_config_column


class PrecisionCalculator:
    def __init__(self):
        # Initialize an empty DataFrame with the required columns
        self.results = pd.DataFrame(columns=["Category", "Precision", "Recall", "F1-score", "TP", "FP", "FN"])

    def add_confusion_matrix(self, type_name, tp, fp, fn):
        Precision = tp / (tp + fp) if (tp + fp) > 0 else 0
        Recall = tp / (tp + fn) if (tp + fn) > 0 else 0
        F1Score = 2 * Precision * Recall / (Precision + Recall) if (Precision + Recall) > 0 else 0

        # Create a new row as a DataFrame and append it to self.results
        new_row = pd.DataFrame([{
            "Category": type_name,
            "Precision": Precision,
            "Recall": Recall,
            "F1-score": F1Score,
            "TP": tp,
            "FP": fp,
            "FN": fn
        }])

        self.results = pd.concat([self.results, new_row], ignore_index=True)
        return self

    def write_results(self, output_path: str):
        self.results.to_csv(output_path, index=False, sep='\t')
        return self

    def print_confusion_matrix_metrics(self, mode="full"):
        """
        Prints the results. 
        Modes:
        - "full": Prints all metrics (Category, Precision, Recall, F1-score, TP, FP, FN).
        - "summary": Prints only Category, Precision, TP, FP.
        """
        if mode == "full":
            # Full format including all metrics
            header_format = "{:<50} {:<9} {:<9} {:<9} {:<5} {:<5} {:<5}"
            row_format = "{:<50} {:<9.4f} {:<9.4f} {:<9.4f} {:<5d} {:<5d} {:<5d}"
            header = header_format.format("Category", "Precision", "Recall", "F1-score", "TP", "FP", "FN")
        elif mode == "summary":
            # Summary format with fewer metrics
            header_format = "{:<35} {:<9} {:<5} {:<5}"
            row_format = "{:<35} {:<9.4f} {:<5d} {:<5d}"
            header = header_format.format("Category", "Precision", "TP", "FP")
        else:
            raise ValueError("Invalid mode. Use 'full' or 'summary'.")

        # Store the results for printing
        results_lines = [header]

        # Iterate through DataFrame rows
        for _, result in self.results.iterrows():
            if mode == "full":
                results_lines.append(
                    row_format.format(
                        result["Category"],
                        float(result["Precision"]),
                        float(result["Recall"]),
                        float(result["F1-score"]),
                        int(result["TP"]),
                        int(result["FP"]),
                        int(result["FN"])
                    )
                )
            elif mode == "summary":
                results_lines.append(
                    row_format.format(
                        result["Category"],
                        float(result["Precision"]),
                        int(result["TP"]),
                        int(result["FP"])
                    )
                )

        # Combine results and print them
        results_output = "\n".join(results_lines)
        logging.info("\n\n" + results_output + "\n")
        return self

    def clear_results(self):
        # Reset the DataFrame to empty with the same columns
        self.results = pd.DataFrame(columns=["Category", "Precision", "Recall", "F1-score", "TP", "FP", "FN"])
        return self

    def sort_results(self, column: str = "Category", reverse: bool = False):
        # Sort the DataFrame
        ascending = not reverse
        self.results = self.results.sort_values(by=column, ascending=ascending).reset_index(drop=True)
        return self

    def output_picture(self, output_path: str, format: list = ["sample", "purity", "software", "variant_type"]):
        """
        Generate a single combined visualization plot for Precision, Recall, and F1-score across different samples and variant types.
        
        Args:
            output_path (str): Directory path to save the plots
            format (list): List of column names to split from Category column
        """
        # Create output directory if not exists
        plot_dir = f"{output_path}/metrics_plot"
        os.makedirs(plot_dir, exist_ok=True)
        
        # Split Category column and convert purity to numeric
        self.results[format] = self.results["Category"].str.split(",", expand=True).apply(lambda x: x.str.strip())
        self.results["purity"] = pd.to_numeric(self.results["purity"], errors='coerce')
        self.results.query("purity != 1", inplace=True)
        print(self.results["software"].unique())
        
        # filter_softwares = ['ClairS_TO_v0_3_0_pileup_nonsomatic', 'Longphase_TO_v0_0_1', 'ClairS_TO_v0_3_0', 'DeepSomatic_TO_v1_8_0']
        #         ['DeepSomatic_TO_v1_8_0' 'ClairS_TO_v0_3_0'
        #  'ClairS_TO_v0_3_0_pileup_nonsomatic' 'Longphase_TO_v0_0_1'
        #  'Longphase_TO_v0_0_1_2' 'Longphase_TO_v0_0_1_deepPON'
        #  'Longphase_TO_v0_0_1_deepVCF']
        filter_softwares = ['ClairS_TO_ssrs_v0_3_0_pileup']
        filter_softwares.extend(['ClairS_TO_ssrs_v0_3_0_pileup_nonsomatic'])
        filter_softwares.extend(['ClairS_TO_ssrs_v0_3_0'])
        filter_softwares.extend(['Longphase_TO_ssrs_v0_0_1'])
        filter_softwares.extend(['Longphase_TO_ss_v0_0_1'])
        # filter_softwares.extend(['Longphase_TO_deep_v0_0_1'])
        # filter_softwares.extend(['Longphase_TO_deep_pileup_v0_0_1'])
        filter_softwares.extend(['ClairS_TO_ss_v0_3_0'])
        filter_softwares.extend(['ClairS_TO_ss_v0_3_0_pileup'])
        filter_softwares.extend(['ClairS_TO_ss_v0_3_0_pileup_nonsomatic'])
        filter_softwares.extend(['DeepSomatic_TO_v1_8_0'])
        # filter_softwares = []

        # filter_softwares.extend(['Longphase_TO_v0_0_1_deepVCF', 'Longphase_TO_v0_0_1_deepPON'])
        # # filter_softwares.extend(['Longphase_TO_v0_0_1_deepPON'])
        # filter_softwares.extend(['Longphase_TO_v0_0_1_deepvcf'])
        # filter_softwares.extend(['Longphase_TO_v0_0_1_3_second75'])
        # filter_softwares.extend(['Longphase_TO_v0_0_1_second75'])
        # filter_softwares.extend(['Longphase_TO_v0_0_1_cluster_nhom_nremove'])
        # filter_softwares.extend(['Longphase_TO_v0_0_1_cluster_nhom'])
        # # filter_softwares.extend(['Longphase_TO_v0_0_1_diss3-add8'])
        # # filter_softwares.extend(['Longphase_TO_v0_0_1'])
        # filter_softwares.extend(['ClairS_TO_v0_3_0', 'DeepSomatic_TO_v1_8_0'])
        # filter_softwares.extend(['Longphase_TO_v0_0_1_2', 'Longphase_TO_v0_0_1_1'])
        # # filter_softwares.extend(['Longphase_TO_v0_0_1_curr'])
        # # filter_softwares.extend(['Longphase_TO_v0_0_1_1'])
        # filter_softwares.extend(['Longphase_TO_v0_0_1', 'Longphase_TO_v0_0_1_1'])
        self.results = self.results.query("software not in @filter_softwares")
        # self.results["software"] = self.results["software"].replace("Longphase_TO_v0_0_1", "Longphase_TO_v0_0_1_3")
        # self.results["Category"] = self.results["Category"].str.replace("Longphase_TO_v0_0_1", "Longphase_TO_v0_0_1_3", regex=False)
        # self.results["software"] = self.results["software"].replace("Longphase_TO_v0_0_1_2", "Longphase_TO_v0_0_1")
        # self.results["Category"] = self.results["Category"].str.replace("Longphase_TO_v0_0_1_2", "Longphase_TO_v0_0_1", regex=False)
        # self.results["software"] = self.results["software"].replace("Longphase_TO_v0_0_1_curr", "Longphase_TO_v0_0_1")
        # self.results["Category"] = self.results["Category"].str.replace("Longphase_TO_v0_0_1_curr", "Longphase_TO_v0_0_1", regex=False)
        metric_list = ["F1-score", "Precision", "Recall"]
        
        # Add faceting plot generation
        for variant_type in self.results["variant_type"].unique():
            variant_data = self.results[self.results["variant_type"] == variant_type]
            samples = sorted(variant_data["sample"].unique(), reverse=False)
            softwares_reverse = sorted(variant_data["software"].unique(),
                            key=lambda x: software_rank.get(x, software_rank['default']),
                            reverse=True)
            softwares_ascending = sorted(variant_data["software"].unique(),
                            key=lambda x: software_rank.get(x, software_rank['default']),
                            reverse=False)
            
            # Calculate grid dimensions
            n_rows = len(metric_list)
            n_cols = len(samples)
            
            # Create figure with subplots
            # fig, axes = plt.subplots(n_rows, n_cols, figsize=(3*n_cols, 15*n_rows), sharey='row')
            fig, axes = plt.subplots(n_rows, n_cols, figsize=(3*n_cols, 4.5*n_rows), sharey='row')
            fig.suptitle(f"Metrics by Tumor Purity Across Samples – {variant_type}", fontsize=24)
            marker_list = {'default': markers.pop(5)}  # Remove and use the marker at index 5
            for software in softwares_reverse:
                marker = markers.pop(random.randint(0, len(markers) - 1))  # Randomly select a marker
                marker_list[software] = marker
            # Plot each metric and sample combination
            for row_idx, metric in enumerate(metric_list):
                for col_idx, sample in enumerate(samples):
                    ax = axes[row_idx, col_idx] if n_rows > 1 else axes[col_idx]
                    sample_data = variant_data[variant_data["sample"] == sample]
                    
                    # Plot data for each software
                    for software in softwares_reverse:
                        if software not in software_colors:
                            software_colors[software] = get_random_color()
                        software_data = sample_data[sample_data["software"] == software].sort_values("purity")
                        color = software_colors.get(software, software_colors['default'])
                        line_style = line_styles.get(software, line_styles['default'])
                        marker = marker_list.get(software, marker_list['default'])
                        
                        ax.plot(software_data["purity"], 
                            software_data[metric],
                            marker=marker,
                            color=color,
                            linewidth=2,
                            markersize=6,
                            label=software)
                    
                    # Customize subplot
                    if row_idx == 0:
                        ax.set_title(f"{sample}", fontsize=16)
                    if row_idx != n_rows - 1:
                        ax.tick_params(labelbottom=False)  # 隱藏非最後 row 的 x tick labels
                        ax.tick_params(axis='x', length=0)      # 不要刻度的小短線
                    if row_idx == n_rows - 1:
                        ax.set_xticks(software_data["purity"].unique())
                        ax.set_xlabel("Purity", fontsize=16)
                    if col_idx == 0:
                        ax.set_ylabel(f"{metric}", fontsize=20)
                    if col_idx != 0:
                        ax.tick_params(labelleft=False)
                        # # ax.set_yticks([])                       # 不要刻度值的位置
                        # ax.set_yticklabels([])                  # 不要刻度的文字
                        # ax.set_ylabel("")                       # 不要 y 軸 label
                        # # ax.spines["left"].set_visible(False)    # 不要左邊的軸線（框線）
                        ax.tick_params(axis='y', length=0)      # 不要刻度的小短線
                        # # ax.yaxis.grid(False)                    # 不要 y 軸的網格線
                    # if idx != len(samples) - 1:
                    #     ax.spines["right"].set_visible(False)    # 不要右邊的軸線（框線）
                    ax.set_ylim(0, 1)
                    ax.set_yticks(np.arange(0, 1.01, 0.1))
                    ax.grid(True, linestyle='--', alpha=0.7)
                    
            
            # Add shared legend
            handles = [plt.Line2D([0], [0],
                                  marker=marker_list.get(sw, marker_list['default']),
                                  color=software_colors.get(sw, software_colors['default']),
                                  label=sw,
                                  linewidth=2)
                       for sw in softwares_ascending]
            fig.legend(handles=handles, title="Software", loc='center left', bbox_to_anchor=(0.95, 0.5),
                        fontsize=14,             # Increase label font size (e.g., 14)
                        title_fontsize=16,       # Increase title font size (e.g., 16)
                        markerscale=2)         # Make markers 1.5x larger in the legend
                            
            # Adjust layout and save
            plt.tight_layout(rect=[0, 0, 0.95, 0.94])
            plt.savefig(f"{plot_dir}/{variant_type}_metrics_combined.png", 
                    dpi=300, 
                    bbox_inches='tight')
            plt.close()
            print(f"Saved {plot_dir}/{variant_type}_metrics_combined.png")
            
        return self

    def _create_legend_handles(self, samples: list, softwares: list) -> tuple:
        """Create legend handles for samples and software."""
        sample_handles = []
        software_handles = []
        
        for idx, sample in enumerate(samples):
            sample_handles.append(plt.Line2D([0], [0], 
                                           marker=markers[idx % len(markers)],
                                           color='black',
                                           label=sample,
                                           markerfacecolor='black',
                                           markersize=8,
                                           linestyle='None'))
        
        for software in softwares:
            color = software_colors.get(software, software_colors['default'])
            line_style = line_styles.get(software, line_styles['default'])
            software_handles.append(plt.Line2D([0], [0],
                                             color=color,
                                             linestyle=line_style,
                                             label=software,
                                             linewidth=2))
        
        return sample_handles, software_handles
    
    def _create_x_axis_labels(self, variant_data: pd.DataFrame, samples: list) -> tuple:
        """Create x-axis labels and positions with two levels: purity on top and sample in the middle.
        
        Returns:
            tuple: (x_labels, x_positions, mapping) where mapping is (sample, purity) -> x_position
        """
        x_labels = []
        x_positions = []
        mapping = {}
        current_pos = 0
        
        for sample in samples:
            sample_data = variant_data[variant_data["sample"] == sample]
            purities = sorted(sample_data["purity"].unique())
            n = len(purities)
            mid_index = n // 2
            
            for i, purity in enumerate(purities):
                # Create mapping
                mapping[(sample, purity)] = current_pos
                
                # Only show sample label in the middle position
                if i == mid_index:
                    label = f"{purity}\n{sample}"
                else:
                    label = f"{purity}\n"
                x_labels.append(label)
                x_positions.append(current_pos)
                current_pos += 1
                
        return x_labels, x_positions, mapping
    
    def _customize_and_save_plot(self, plot_dir: str, filename: str, title: str,
                               sample_handles: list = None, software_handles: list = None,
                               x_labels: list = None, x_positions: list = None):
        """Customize and save the plot."""
        plt.xlabel("Purity", fontsize=12)
        plt.ylabel("F1 Score", fontsize=12)
        plt.title(title, fontsize=14)
        plt.grid(True, linestyle='--', alpha=0.7)
        
        if sample_handles and software_handles:
            legend1 = plt.legend(handles=sample_handles, 
                               title="Sample",
                               bbox_to_anchor=(1.15, 1),
                               loc='upper left')
            plt.gca().add_artist(legend1)
            plt.legend(handles=software_handles,
                      title="Software",
                      bbox_to_anchor=(1.15, 0.5),
                      loc='center left')
        else:
            plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        
        if x_labels and x_positions:
            plt.xticks(x_positions, x_labels, rotation=45, ha='right')
        
        plt.ylim(0, 1)
        plt.tight_layout()
        
        plt.savefig(f"{plot_dir}/{filename}", 
                   dpi=300, 
                   bbox_inches='tight')
        plt.close()

    def output_radar_picture(self, output_path: str, format: list = ["sample", "purity", "software", "variant_type"]):
        if not os.path.exists(f"{output_path}/metrics_plot"):
            os.makedirs(f"{output_path}/metrics_plot")
        self.results[format] = self.results["Category"].str.split(",", expand=True).apply(lambda x: x.str.strip())
        variant_types = self.results["variant_type"].unique()  # Extract unique variant types
        
        pbar = tqdm(total=(len(variant_types)*3*5*2+1), desc="Processing")
        for variant_type in variant_types:
            for picture_name in ["Precision", "Recall", "F1-score"]:
                # Get all purity values
                purities = self.results["purity"].unique()
                
                # Create a radar chart for each purity value with progress bar
                for purity in purities:
                    pbar.set_description(f"Executing: {variant_type} {picture_name} {purity}")
                    pbar.update(1)
                    # 將單獨雷達圖尺寸調整到較大的尺寸，並設置更高dpi以呈現更多細節
                    plt.figure(figsize=(12, 8))
                    # plt.figure(figsize=(18, 12))
                    
                    # 過濾當前變異類型和純度的數據
                    filtered_data = self.results[(self.results["variant_type"] == variant_type) & 
                                                (self.results["purity"] == purity)]
                    
                    # 獲取唯一的樣本和軟件
                    samples = filtered_data["sample"].unique()
                    softwares = filtered_data["software"].unique()
                    
                    # 如果沒有數據，跳過此圖
                    if len(samples) == 0 or len(softwares) == 0:
                        plt.close()
                        continue
                    
                    # 設置雷達圖
                    angles = np.linspace(0, 2*np.pi, len(samples), endpoint=False).tolist()
                    angles += angles[:1]  # 閉合圖形
                    
                    # 設置圖表
                    ax = plt.subplot(111, polar=True)
                    
                    # 添加每個樣本的標籤（角標）
                    plt.xticks(angles[:-1], samples, fontsize=12)
                    
                    # 設置y軸限制
                    ax.set_ylim(0, 1)
                    
                    # 為每個軟件繪製一條線
                    for software in softwares:
                        software_data = filtered_data[filtered_data["software"] == software].copy()
                        
                        # 創建一個字典，將樣本映射到指標值
                        sample_to_value = dict(zip(software_data["sample"], software_data[picture_name]))
                        
                        # 按照樣本順序獲取值
                        values = [sample_to_value.get(sample, 0) for sample in samples]
                        values += values[:1]  # 閉合圖形
                        
                        # 繪製線條
                        ax.plot(angles, values, linewidth=2, label=software)
                        ax.fill(angles, values, alpha=0.25)
                    
                    # 添加圖例和標題
                    plt.legend(loc='upper right', bbox_to_anchor=(0.1, 0.1), fontsize=12)
                    plt.title(f"{picture_name} Radar Chart for {variant_type} (Purity: {purity})", fontsize=16)
                    
                    # 保存單獨圖表，並設置高dpi
                    plt.savefig(f"{output_path}/metrics_plot_radar/{variant_type}_{picture_name}_{purity}.png", dpi=300)
                    # logging.info(f"Saved {output_path}/metrics_plot/{variant_type}_{picture_name}_{purity}.png")
                    plt.close()
                
            # 合併所有圖表為一個大圖表
            # 按照指標（Precision、Recall、F1-score）和純度（purity）排列
            metrics = ["Precision", "Recall", "F1-score"]
            purities = sorted(self.results["purity"].unique())
            
            # 創建合併圖時，將每個子圖設置為較大的尺寸
            fig, axes = plt.subplots(len(metrics), len(purities), figsize=(len(purities)*6, len(metrics)*5))
            
            # 設置圖表標題
            fig.suptitle(f"Combined Radar Charts for {variant_type}", fontsize=18)
            # 添加行標籤（指標）
            for i, metric in enumerate(metrics):
                axes[i, 0].set_ylabel(metric, fontsize=16, rotation=90)
            
            # 添加列標籤（純度）
            for j, purity in enumerate(purities):
                axes[0, j].set_title(f"Purity: {purity}", fontsize=16)

            # 填充圖表
            for i, metric in enumerate(metrics):
                for j, purity in enumerate(purities):
                    pbar.set_description(f"Merging: {variant_type} {metric} {purity}")
                    pbar.update(1)
                    img_path = f"{output_path}/metrics_plot/{variant_type}_{metric}_{purity}.png"
                    if os.path.exists(img_path):
                        img = plt.imread(img_path)
                        axes[i, j].imshow(img)
                        axes[i, j].axis('off')
                    else:
                        axes[i, j].text(0.5, 0.5, "No Data", ha='center', va='center', fontsize=16)
                        axes[i, j].axis('off')
            pbar.set_description(f"Writing")
            plt.tight_layout(rect=[0, 0, 1, 0.96])
            plt.savefig(f"{output_path}/{variant_type}_metrics.png", dpi=300, bbox_inches='tight', pad_inches=0)
            plt.close()
            pbar.update(1)
            logging.info(f"Saved combined radar chart: {output_path}/{variant_type}_metrics.png")
            
        return self

def get_result(precision_calculator, name, in_vcf_path=None, in_vcf_df=None, ans_vcf_path=None, ans_bed_path=None, test_chr=None, use_alt=True, check_type="snp", filter_type="PASS", no_filter=False):
    if check_type == "snv":
        check_type = "snp"
    merge_cols = ["CHROM", "POS", "ALT"] if use_alt else ["CHROM", "POS"]

    # 讀取答案
    ans_vcf_df = BioinfoFileReader(ans_vcf_path).reader().check_variant_type().df.query(f"type_simple == '{check_type}'")[merge_cols]

    # 讀取輸入
    vcf_df = None
    if in_vcf_df is None and ".vcf" in in_vcf_path:
        vcf_df = BioinfoFileReader(in_vcf_path).reader().check_variant_type().df
        if no_filter:
            vcf_df = vcf_df.query(f"type_simple == '{check_type}'")
        else:
            vcf_df = vcf_df.query(f"FILTER == '{filter_type}' and type_simple == '{check_type}'")
        vcf_df = vcf_df.drop_duplicates(subset=["CHROM", "POS"], keep='first')
        # vcf_df = vcf_df[merge_cols]
    else:
        use_alt = False
        vcf_df = in_vcf_df
        vcf_df["POS"] = vcf_df["POS"] + 1

    # 特定比對
    if test_chr is not None:
        ans_vcf_df.query(f"CHROM == '{test_chr}'", inplace=True)
        vcf_df.query(f"CHROM == '{test_chr}'", inplace=True)
    # 比對
    vcf_df = ans_tag(vcf_df, ans_vcf_df, use_alt=use_alt)

    # 讀bed
    if ans_bed_path is not None and os.path.exists(ans_bed_path):
        ans_bed_df = BioinfoFileReader(ans_bed_path).reader().df
        vcf_df = Diff_bed_vcf(ans_bed_df, vcf_df, 3).find_range_muti_chrs("HIGH").vcf_df.query("HIGH == 1")

    TP = len(vcf_df.query("ans == 'tp'"))
    FP = len(vcf_df.query("ans == 'fp'"))
    FN = len(vcf_df.query("ans == 'fn'"))
    precision_calculator.add_confusion_matrix(name, TP, FP, FN)

    return vcf_df


class PhaseResultPrinter:
    def __init__(self):
        self.row_format = "{:<20} {:<10d} {:<13.2f} {:<10d} {:<10.2f} {:<20.2f} {:<12d} {:<13d} {:<10d}"
        self.results_lines = []
        run_root_path = get_config_column("run_root_path")
        self.ans_vcf_path = f"{run_root_path}/HG002_GRCh38_1_22_v4.2.1_benchmark_hifiasm_v11_phasetransfer.dirty.vcf"
        self.whatshap_path = f"{run_root_path}/miniconda3/envs/whatshap/bin/whatshap"

    def add_results_lines_header(self):
        header_format = "{:<20} {:<10} {:<13} {:<10} {:<10} {:<20} {:<12} {:<13} {:<10}"
        header = header_format.format("Category", "Phased_SNP", "Phased_SNP(%)", "SW", "SW(%)", "Hamming_distance(%)", "No._of_Block", "Block_N50(bp)", "Block_Sum")
        self.results_lines.insert(0, header)

    def make_command(self,
                     diff_path: str,
                     ans_vcf_path: str = None,
                     whatshap_path: str = None):
        if ans_vcf_path is None:
            ans_vcf_path = self.ans_vcf_path
        if whatshap_path is None:
            whatshap_path = self.whatshap_path

        # 構造各個命令字串
        make_compare_command = f"{whatshap_path} compare {ans_vcf_path} {diff_path}.vcf --switch-error-bed {diff_path}.bed > {diff_path}.compare"
        make_block_command = f"""
grep -e "0|1" -e "1|0" {diff_path}.vcf |
grep -v "#" |
awk '{{
    split($9, tags, ":")
    split($10, values, ":")
    for (i=1; i<=length(tags); i++) {{
        if (tags[i] == "GT") {{
            gt_value = values[i]
        }}
        if (tags[i] == "PS") {{
            ps_value = values[i]
        }}
    }}
    print $1 "\\t" $2 "\\t" gt_value "\\t" ps_value
}}' |
sort -k1,1 -k4,4n -k2,2n |
awk  '{{
    if( $4 != upPS ){{
        if(NR>1 && count > 1){{
            print upChr "\\t" strPos "\\t" upPos - strPos "\\t" count
        }}

        strPos = $2
        count=0
    }}
    upPS = $4
    upPos = $2
    upChr = $1
    count++
}}
END{{
    print upChr "\\t" strPos "\\t" upPos - strPos "\\t" count
}}
' | sort -k3rn >{diff_path}.block
"""
        phased_count_command = f"awk '!/^#/ {{phased += gsub(/\\|/, \"\")}} END {{print phased}}' {diff_path}.vcf"
        unphased_count_command = f"awk '!/^#/ && length($5) == 1 && length($6) == 1 {{unphased += gsub(/0\\/1:/, \"\")}} END {{print unphased}}' {diff_path}.vcf"
        switch_error_command = f"""
grep -e "Chromosome" -e "phased pairs of variants assessed:" -e "switch errors:" "{diff_path}.compare" |
awk 'NR%5<=3 && NR%5>=1' |
awk '{{
if(NR%3==1){{printf $2"\\t"$3}} 
if(NR%3==2){{printf "\\t"$6}} 
if(NR%3==0){{printf "\\t"$3"\\n"}}}}' |
awk '{{total+=$3;error+=$4;}}END{{print total, error, error/total}}'
"""
        hamming_command = f"""
grep -e "Chromosome" -e "--> covered variants:" -e "Block-wise Hamming distance:" {diff_path}.compare |
awk 'NR%5==1 || NR%5==4 || NR%5==0' |
awk '{{if(NR%3==1){{printf $3"\\t"}}if(NR%3==2){{printf $4"\\t"}}if(NR%3==0){{printf $4"\\n"}}}}' |
awk '{{variant+=$2;Hamming+=$3}}END{{print variant, Hamming, Hamming/variant}}'
"""
        block_N50_command = f"""
awk 'ARGIND==1 {{ sum+=$3 }} 
     ARGIND==2 {{ count+=$3; 
                if(count >= sum/2) {{print $3; exit}} }}' {diff_path}.block {diff_path}.block
"""
        block_count_command = f"""
wc -l {diff_path}.block | awk '{{print $1}}'
"""
        block_sum_command = f"""
awk '{{x=x+$3}} END{{print x}}' {diff_path}.block
"""

        # 將所有命令和相對應需要存儲的結果鍵放到一個 list 中
        self.command_and_args_list = [
            ["make_compare", make_compare_command, []],
            ["make_block", make_block_command, []],
            ["phased_count", phased_count_command, ['phased']],
            ["unphased_count", unphased_count_command, ['unphased']],
            ["switch_error", switch_error_command, ['sw_total', 'sw_error', 'sw_error_rate']],
            ["hamming", hamming_command, ['variant', 'hamming', 'hamming_variant_ratio']],
            ["block_N50", block_N50_command, ['n50_value']],
            ["block_count", block_count_command, ['block_count']],
            ["block_sum", block_sum_command, ['block_sum']]
        ]

    def execute_commands(self,
                         result_dict: dict,
                         show_progress: bool = True,
                         log_file_path: str = "execution_log.txt",
                         check_return: bool = True):

        if show_progress:
            pbar = tqdm(total=len(self.command_and_args_list), desc="Processing Commands")
        else:
            pbar = None

        # 執行命令並寫入日誌檔案
        with open(log_file_path, "w") as log_file:
            for command_and_args in self.command_and_args_list:
                command_name = command_and_args[0]
                command_str = command_and_args[1]
                args = command_and_args[2]
                if show_progress:
                    pbar.set_description(f"Executing: {command_name}")
                    pbar.update(1)
                log_file.write(f"{command_str}\n")

                try:
                    process = subprocess.run(command_str, shell=True, check=check_return, capture_output=True, text=True)
                except subprocess.CalledProcessError as e:
                    log_file.write(f"Command '{command_str}' failed with return code {e.returncode}\n")
                    logging.error(f"Command '{command_str}' failed with return code {e.returncode}")
                    continue

                if args:
                    output_values = process.stdout.strip().split()
                    result_dict.update({key: value for key, value in zip(args, output_values)})

    def process_result(self, result_dict):
        # 計算 phase_rate
        phased = int(result_dict["phased"])
        unphased = int(result_dict["unphased"])
        phase_rate = float(phased / (phased + unphased)) if (phased + unphased) > 0 else 0
        result_dict["phase_rate"] = phase_rate

        # 格式化結果輸出
        self.results_lines.append(
            self.row_format.format(
                result_dict["category"],
                int(result_dict["phased"]),
                float(result_dict["phase_rate"]) * 100,
                int(result_dict["sw_error"]),
                float(result_dict["sw_error_rate"]) * 100,
                float(result_dict["hamming_variant_ratio"]) * 100,
                int(result_dict["block_count"]),
                int(result_dict["n50_value"]),
                int(result_dict["block_sum"])
            )
        )

    def get_result_dict(self,
                        category: str,
                        diff_path: str,
                        ans_vcf_path: str = None,
                        whatshap_path: str = None,
                        show_progress: bool = True,
                        log_file_path: str = "execution_log.txt",
                        check_return: bool = True):
        if ans_vcf_path is None:
            ans_vcf_path = self.ans_vcf_path
        if whatshap_path is None:
            whatshap_path = self.whatshap_path

        # 初始化結果字典
        result_dict = {
            "category": category,
            'phased': 0,
            'unphased': 0,
            # 'phase_rate': 0,
            'sw_total': 0,
            'sw_error': 0,
            'sw_error_rate': 0,
            'variant': 0,
            'hamming': 0,
            'hamming_variant_ratio': 0,
            'n50_value': 0,
            'block_count': 0,
            'block_sum': 0
        }

        self.make_command(diff_path, ans_vcf_path, whatshap_path)
        self.execute_commands(result_dict, show_progress, log_file_path, check_return)
        # self.process_results(result_dict)
        return result_dict

    def get_result_lines(self):
        """返回結果字典"""
        return self.results_lines

    def print_results(self, sort_by_category: bool = True):
        """直接打印結果表格"""
        if sort_by_category:
            self.results_lines.sort(key=lambda x: x.split()[0])
        self.add_results_lines_header()
        results_output = "\n".join(self.results_lines)
        logging.info("\n" + results_output + "\n")
