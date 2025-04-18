from handler.utils import *
from handler.config_parser import get_config_column, get_my_address, get_boos_address, get_boos_default_path
from handler.reader import *
from handler.metrics import *
from handler.command_runner import *
from handler.purity_prediction import *
from purity_ascat.test import *
from purity_ascat.test2 import *
from purity_ascat.test3 import *
from handler.image_splitter import run_image_splitter
def test3(plot_type='boxplot'):
# def test3(plot_type='violinplot'):
    """
    绘制LOH和non-LOH区域的变异等位基因频率分布图
    
    参数:
        plot_type: 图表类型，可选 'boxplot' 或 'violinplot'，默认为 'boxplot'
    """
    data_root_path = get_config_column("data_root_path")
    run_root_path = get_config_column("run_root_path")
    output_folder = get_config_column("output_folder")
    longphase_vcf = f'{data_root_path}/HCC1395/ONT/subsample/t50_n00/ClairS_TO_v0_3_0/snv.vcf'
    one_loh_bed = f"{run_root_path}/{output_folder}/HCC1395_ONT/t50_n00/output_LOH.bed"
    # one_small_cnv_bed = f"{run_root_path}/HCC1395/ONT/subsample/t50_n00/ClairS_TO_v0_3_0/small_cnv.bed"
    bed_loh_df = BioinfoFileReader(one_loh_bed).reader().df
    vcf_df = BioinfoFileReader(longphase_vcf).reader().cut_format_sample(['AF']).df
    vcf_df['AF'] = vcf_df['AF'].astype(float).round(2)
    vcf_df = vcf_df.query("CHROM == 'chr2'")
    vcf_df = Diff_bed_vcf(bed_loh_df, vcf_df, 3).find_range_muti_chrs("LOH").vcf_df
    af_1_df = vcf_df[vcf_df['LOH'] == 1]
    af_0_df = vcf_df[vcf_df['LOH'] == 0]
    
    # 创建输出目录
    output_dir = f"{run_root_path}/test-plot"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # LOH 图表
    plt.figure(figsize=(6, 6))
    
    if plot_type == 'boxplot':
        # 箱线图
        plt.boxplot(
            af_1_df['AF'],
            vert=True,
            labels=['LOH'],
            widths=0.5,
            patch_artist=True,  # 启用填充功能
            boxprops=dict(facecolor='palegreen', color='green'),  # 设置内部颜色和边框颜色
            medianprops=dict(color='darkgreen'),  # 设置中位数线条颜色
            whiskerprops=dict(color='green'),  # 设置须线颜色
            capprops=dict(color='green'),  # 设置须端颜色
            flierprops=dict(markerfacecolor='darkgreen', marker='o')  # 设置异常值样式
        )
    else:
        # 小提琴图
        violin_parts = plt.violinplot(
            af_1_df['AF'],
            vert=True,
            showmeans=True,
            showmedians=True,
            showextrema=True
        )
        
        # 设置小提琴图的颜色
        for pc in violin_parts['bodies']:
            pc.set_facecolor('palegreen')
            pc.set_edgecolor('green')
            pc.set_alpha(0.7)
        
        # 设置统计线的颜色
        violin_parts['cmeans'].set_color('darkgreen')
        violin_parts['cmedians'].set_color('darkgreen')
        violin_parts['cmaxes'].set_color('green')
        violin_parts['cmins'].set_color('green')
        violin_parts['cbars'].set_color('green')
        
        plt.xticks([1], ['LOH'])
    plt.title(f'LOH {plot_type}')
    output_file = f"{output_dir}/af_{plot_type}_loh.png"
    
    plt.ylabel('AF')
    
    # 添加统计信息
    # mean_af = af_1_df['AF'].mean()
    # median_af = af_1_df['AF'].median()
    # plt.text(0.05, 0.95, f'Mean: {mean_af:.3f}\nMedian: {median_af:.3f}', 
    #          transform=plt.gca().transAxes, verticalalignment='top',
    #          bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()
    
    # non-LOH 图表
    plt.figure(figsize=(6, 6))
    
    if plot_type == 'boxplot':
        # 箱线图
        plt.boxplot(
            af_0_df['AF'],
            vert=True,
            labels=['non-LOH'],
            widths=0.5,
            patch_artist=True,  # 启用填充功能
            boxprops=dict(facecolor='skyblue', color='blue'),  # 设置内部颜色和边框颜色
            medianprops=dict(color='darkblue'),  # 设置中位数线条颜色
            whiskerprops=dict(color='blue'),  # 设置须线颜色
            capprops=dict(color='blue'),  # 设置须端颜色
            flierprops=dict(markerfacecolor='darkblue', marker='x')  # 设置异常值样式
        )
    else:
        # 小提琴图
        violin_parts = plt.violinplot(
            af_0_df['AF'],
            vert=True,
            showmeans=True,
            showmedians=True,
            showextrema=True
        )
        
        # 设置小提琴图的颜色
        for pc in violin_parts['bodies']:
            pc.set_facecolor('skyblue')
            pc.set_edgecolor('blue')
            pc.set_alpha(0.7)
        
        # 设置统计线的颜色
        violin_parts['cmeans'].set_color('darkblue')
        violin_parts['cmedians'].set_color('darkblue')
        violin_parts['cmaxes'].set_color('blue')
        violin_parts['cmins'].set_color('blue')
        violin_parts['cbars'].set_color('blue')
        
        plt.xticks([1], ['non-LOH'])
    plt.title(f'non-LOH {plot_type}')
    output_file = f"{output_dir}/af_{plot_type}_non_loh.png"
    
    plt.ylabel('AF')
    
    # 添加统计信息
    # mean_af = af_0_df['AF'].mean()
    # median_af = af_0_df['AF'].median()
    # plt.text(0.05, 0.95, f'Mean: {mean_af:.3f}\nMedian: {median_af:.3f}', 
    #          transform=plt.gca().transAxes, verticalalignment='top',
    #          bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()
    
    print(f"已生成 {plot_type} 图表，保存在 {output_dir} 目录下")


def test6():
    data_root_path = get_config_column("data_root_path")
    path = f"{data_root_path}/orthogonal-tools-benchmark/vcfs/"
    vcf_path = path + "H1437_orthogonal-tools-benchmark_somatic-only.vcf"
    bed_path = path.replace("vcfs", "beds") + \
        "H1437_orthogonal-tools-benchmark.bed"
    vcf_df = BioinfoFileReader(vcf_path).reader().check_variant_type().df.query("type_simple == 'snp'").drop(columns=["type_simple", "type_complex"]).drop(
        columns=["TUMOR", "2:NORMAL", "2:TUMOR", "2:ClairS", "3:ClairS", "3:2:ClairS"]).rename(columns={"NORMAL": "SAMPLE"})
    bed_df = BioinfoFileReader(bed_path).reader().df
    tag_df = Diff_bed_vcf(bed_df, vcf_df, 3).find_range_muti_chrs(
        "HIGH").vcf_df.query("HIGH == 1")
    tag_df.to_csv(f"{path}H1437.body", sep="\t", index=False, header=False)
    subprocess.run(
        f"grep '^#' {vcf_path} > {path}H1437.vcf", shell=True, check=True)
    subprocess.run(
        f"cat {path}H1437.body >> {path}H1437.vcf", shell=True, check=True)
    subprocess.run(
        f"sed 's#\\./\\.:#0|1:#g' {path}H1437.vcf > {path}H1437_phased.vcf", shell=True, check=True)
    # subprocess.run(r"sed 's#\./\.:#0|1:#g' {path}H1437.vcf > {path}H1437_phased.vcf", shell=True, check=True)

def merge_images():
    """
    Download and merge images based on their matrix position from PMC
    """
    import requests
    from PIL import Image
    import io
    import os
    import re

    # Create directory for downloaded images
    run_root_path = get_config_column("run_root_path")
    output_plot_folder = get_config_column("output_plot_folder")
    output_dir = f"{run_root_path}/{output_plot_folder}/downloaded_images"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Base URL and parameters
    base_url = "/corecgi/tileshop/tileshop.fcgi"
    params = {
        "p": "PMC3",
        "id": "420268",
        "s": "128"
    }
    
    # Dictionary to store images by position
    images = {}
    max_row = 0
    max_col = 0
    
    # Store maximum dimensions for each row and column
    row_heights = {}  # Maximum height for each row
    col_widths = {}   # Maximum width for each column
    
    # Download all images and calculate dimensions
    for r in range(6):  # 6 rows
        for c in range(5):  # 5 columns
            params["r"] = r + 1
            params["c"] = c + 1
            
            # Construct full URL
            url = f"https://www.ncbi.nlm.nih.gov{base_url}"
            
            try:
                response = requests.get(url, params=params)
                if response.status_code == 200:
                    # Save image to memory
                    img = Image.open(io.BytesIO(response.content))
                    images[(r, c)] = img
                    
                    # Update maximum dimensions for rows and columns
                    width, height = img.size
                    row_heights[r] = max(row_heights.get(r, 0), height)
                    col_widths[c] = max(col_widths.get(c, 0), width)
                    
                    # Update maximum row and column indices
                    max_row = max(max_row, r)
                    max_col = max(max_col, c)
                    
                    # Save individual image
                    img.save(f"{output_dir}/tile_{r}_{c}.png")
            except Exception as e:
                print(f"Error downloading image at position ({r}, {c}): {e}")
    
    # Calculate total dimensions
    total_width = sum(col_widths.values())
    total_height = sum(row_heights.values())
    
    # Create new image with calculated dimensions
    merged_image = Image.new('RGB', (total_width, total_height))
    
    # Calculate cumulative dimensions for positioning
    cum_heights = {-1: 0}  # Start with 0 at index -1
    cum_widths = {-1: 0}   # Start with 0 at index -1
    
    for i in range(max_row + 1):
        cum_heights[i] = cum_heights[i-1] + row_heights.get(i, 0)
    for j in range(max_col + 1):
        cum_widths[j] = cum_widths[j-1] + col_widths.get(j, 0)
    
    # Paste all images in correct positions
    for (r, c), img in images.items():
        x = cum_widths[c-1]  # Start at the cumulative width of previous columns
        y = cum_heights[r-1] # Start at the cumulative height of previous rows
        merged_image.paste(img, (x, y))
    
    # Save final merged image
    merged_image.save(f"{output_dir}/merged_image.png")
    print(f"Merged image saved to {output_dir}/merged_image.png")

def gui_slip_image():
    """
    Interactive image editor for splitting images with draggable lines
    """
    # Get config and setup paths
    run_root_path = get_config_column("run_root_path")
    output_plot_folder = get_config_column("output_plot_folder")
    
    # input_dir = f"{run_root_path}/{output_plot_folder}/split_images"
    # output_dir = f"{run_root_path}/{output_plot_folder}/cropped_images"
    input_dir = f"{run_root_path}/{output_plot_folder}/transformed_images"
    output_dir = f"{run_root_path}/{output_plot_folder}/transformed_images2"

    # Run the image splitter application
    run_image_splitter(input_dir, output_dir)

def process_image_transform():
    """
    將圖片左右對稱相反過來並向左旋轉90度
    Input: /bip8_disk/zhenyu112/test-plot/cropped_images
    Output: /bip8_disk/zhenyu112/test-plot/transformed_images
    """
    from PIL import Image
    import os

    # 設定輸入和輸出路徑
    run_root_path = get_config_column("run_root_path")
    output_plot_folder = get_config_column("output_plot_folder")
    input_dir = f"{run_root_path}/{output_plot_folder}/cropped_images"
    output_dir = f"{run_root_path}/{output_plot_folder}/transformed_images"

    # 創建輸出目錄
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # 處理所有PNG圖片
    for filename in os.listdir(input_dir):
        if filename.endswith(".png"):
            # 讀取圖片
            image_path = os.path.join(input_dir, filename)
            img = Image.open(image_path)

            # 左右對稱翻轉
            flipped_img = img.transpose(Image.FLIP_LEFT_RIGHT)

            # 向左旋轉90度
            rotated_img = flipped_img.rotate(90, expand=True)

            # 保存處理後的圖片
            output_path = os.path.join(output_dir, f"{filename}")
            rotated_img.save(output_path)
            print(f"已處理並保存: {output_path}")

def slip_image():
    """
    Split image into halves horizontally and 12 parts vertically
    """
    from PIL import Image
    import os

    # Get config and setup paths
    run_root_path = get_config_column("run_root_path")
    output_plot_folder = get_config_column("output_plot_folder")
    input_image_path = f"{run_root_path}/{output_plot_folder}/downloaded_images/merged_image.png"
    output_dir = f"{run_root_path}/{output_plot_folder}/split_images"
    
    # Create output directory if not exists
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Open and load image
    img = Image.open(input_image_path)
    width, height = img.size

    x_part = 12
    y_part = 2
    
    # Calculate dimensions for splitting
    part_width = width // x_part
    part_height = height // y_part
    
    # Split image into sections
    for x in range(x_part):  # 2 horizontal parts
        for y in range(y_part):  # 12 vertical parts
            # Calculate coordinates
            left = x * part_width
            upper = y * part_height
            right = left + part_width
            lower = upper + part_height
            
            # Crop image section
            section = img.crop((left, upper, right, lower))
            
            # Save section
            output_path = f"{output_dir}/section_{x}_{y}.png"
            section.save(output_path)
            
    print(f"Image split into {x_part*y_part} sections and saved to {output_dir}")

def write_image():
    """
    Generate radar charts from metrics file
    """
    run_root_path = get_config_column("run_root_path")
    output_folder = get_config_column("output_folder")
    output_plot_folder = get_config_column("output_plot_folder")
    output_root_path = f"{run_root_path}/{output_plot_folder}"
    input_root_file = f"{run_root_path}/bioinfo_handler/data/metrics.txt"
    
    # Initialize PrecisionCalculator
    calculator = PrecisionCalculator()
    
    # Read metrics file
    if not os.path.exists(input_root_file):
        logging.error(f"Metrics file not found: {input_root_file}")
        return
        
    # Load data into DataFrame
    df = pd.read_csv(input_root_file, sep="\t")
    
    # Convert to desired format and load into calculator
    calculator.results = pd.DataFrame([{
        "Category": row["Category"],
        "Precision": row["Precision"],
        "Recall": row["Recall"],
        "F1-score": row["F1-score"],
        "TP": row["TP"],
        "FP": row["FP"],
        "FN": row["FN"]
    } for _, row in df.iterrows()])
    
    # Split Category column into sample, purity, software, and variant_type
    format_columns = ["sample", "purity", "software", "variant_type"]
    calculator.results[format_columns] = calculator.results["Category"].str.split(",", expand=True)
    calculator.results["purity"] = calculator.results["purity"].astype(float)
    # Convert numeric columns to float
    for col in ["Precision", "Recall", "F1-score"]:
        calculator.results[col] = calculator.results[col].astype(float)
    
    # Convert count columns to int
    for col in ["TP", "FP", "FN"]:
        calculator.results[col] = calculator.results[col].astype(int)
    
    # Generate radar charts
    # logging.info("Generating radar charts...")
    # calculator.output_radar_picture(output_root_path)
    # logging.info(f"Charts have been saved to {output_root_path}/metrics_plot/ and {output_root_path}/")

    # Generate line charts
    logging.info("Generating line charts...")
    calculator.output_picture(output_root_path)
    logging.info(f"Charts have been saved to {output_root_path}/metrics_plot")
    calculator.print_confusion_matrix_metrics(mode="full")

def plot_vaf_distribution(log_file_path, output_dir):
    """
    Read VAF distribution data and generate bar charts
    
    Args:
        log_file_path (str): Input log file path
        output_dir (str): Output image directory
    """
    # Create output directory
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    
    # Set column names
    columns = ['chr', 'start', 'end', 'peak_class', 'threshold', 'num_peaks']
    vaf_columns = [str(i) for i in range(101)]  # 0-100
    columns += vaf_columns
    # Read data
    df = pd.read_csv(log_file_path, sep='\t', header=None, names=columns)
    print(df.head())
    
    # Generate a chart for each chr and interval
    for _, row in df.iterrows():
        chr_name = row['chr']
        start = row['start']
        end = row['end']
        
        # Get VAF data
        vaf_data = row[vaf_columns].astype(float)
        
        # Create figure
        plt.figure(figsize=(12, 6))
        
        # Draw bar chart
        plt.bar(range(101), vaf_data, width=1.0)
        
        # Set title and labels
        plt.title(f'VAF Distribution - Chr{chr_name} ({start}-{end})')
        plt.xlabel('VAF (%)')
        plt.ylabel('Count')
        
        # Set x-axis ticks
        plt.xticks(range(0, 101, 10))
        
        # Add grid
        plt.grid(True, axis='y', linestyle='--', alpha=0.7)
        
        # Adjust layout
        plt.tight_layout()
        
        # Save image
        output_file = os.path.join(output_dir, f'{chr_name}:{start}.png')
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        # print(f"Generated image: {output_file}")
    print(f"Generated {len(df)} images: {output_dir}")

def test_polynomial_regression():
    run_root_path = get_config_column("run_root_path")
    root_dir = f"{run_root_path}/{get_config_column('output_folder')}"
    model_path = f"{root_dir}/test.joblib"
    # 刷新模型
    if os.path.exists(model_path):
        os.remove(model_path)
    degree = 3
    # Get all subdirectories in root_dir
    all_dir = []
    check_dir = []
    
    # Scan all sample directories
    for sample in os.listdir(root_dir):
        sample_dir = os.path.join(root_dir, sample)
        if os.path.isdir(sample_dir) and not sample.endswith('log'):
        # if os.path.isdir(sample_dir) and not sample.endswith('log') and not sample.endswith('Dorado'):
            all_dir.append(sample_dir)

    test_root_dir = f"{run_root_path}/{get_config_column('test_folder')}"
    # Scan all sample directories
    for sample in os.listdir(test_root_dir):
        sample_dir = os.path.join(test_root_dir, sample)
        if os.path.isdir(sample_dir) and not sample.endswith('log'):
        # if os.path.isdir(sample_dir) and not sample.endswith('log') and not sample.endswith('Dorado'):
            all_dir.append(sample_dir)
            check_dir.append(sample_dir)
    # check_dir = all_dir
    train_dir = all_dir

    # Filter out directories containing specified samples
    need_remove = {"HCC1954", "COLO829", "HCC1395"}
    need_remove = {"HCC1954", "COLO829"}
    need_remove = {"HCC1954", "HCC1395"}
    need_remove = {"HCC1954"}
    # need_remove = {}
    train_dir = [d for d in train_dir if not any(sample in d for sample in need_remove)]
    # train_dir.remove(f"{root_dir}/HCC1395_ONT")
    # print(train_dir)

    predictor = PurityPredictor(degree=degree)
    model_info = predictor.train(train_dir, model_path)
    prediction = predictor.predict(check_dir, model_path, plot_model=True, use_function=True)
    # print_grouped_data(prediction)
    if isinstance(prediction, list) and len(prediction) > 0 and isinstance(prediction[0], tuple):
        plot_prediction_vs_actual(prediction, output_path=f"{root_dir}/purity_prediction_results.png")

    # # 訓練模型
    # results, model_info = predict_purity_polynomial_regression(
    #     input_paths=train_dir,
    #     degree=degree,
    #     train_model=True,
    #     model_path=model_path
    # )


    # if len(check_dir) > 0:
    #     # 預測純度
    #     prediction, _ = predict_purity_polynomial_regression(
    #         input_paths=check_dir,
    #         train_model=False,
    #         model_path=model_path,
    #         # plot_model=False
    #         plot_model=True
    #     )
    #     # 印預測結果
    #     print_grouped_data(prediction)

    #     # 如果prediction是一個包含預測結果的列表，則繪製圖表
    #     if isinstance(prediction, list) and len(prediction) > 0 and isinstance(prediction[0], tuple):
    #         plot_prediction_vs_actual(prediction, output_path=f"{root_dir}/purity_prediction_results.png")


def file_for():
    server_name = socket.gethostname()
    data_root_path = get_config_column("data_root_path")
    output_root_path = f"{get_config_column('run_root_path')}/{get_config_column('output_folder')}"
    ans_rank = get_config_column("ans_rank")
    ans_vcf_filter = get_config_column("ans_vcf_filter")
    sequencing_platform_rank = get_config_column("sequencing_platform_rank")
    run_type = get_config_column("run_type")
    run_samples = get_config_column("run_samples")
    run_puritys = get_config_column("run_puritys")
    run_calling_softwares = get_config_column("run_calling_softwares")
    input_calling_softwares = get_config_column("input_calling_softwares")
    plot_vaf_distribution_path = f"{get_config_column('run_root_path')}/{get_config_column('output_plot_folder')}/vaf_distribution"
    tmp_num = get_config_column("tmp")
    # run_type = ["snv"]  # snv indel
    # run_samples = ["H2009", "HCC1954", "HCC1937", "COLO829", "H1437", "HCC1395"]
    # run_samples = ["HCC1395"]
    # run_puritys = [1, 0.8, 0.6, 0.4, 0.2]
    # run_calling_softwares = ["ClairS_TO_v0_3_0", "ClairS_TO_v0_3_0_pileup", "ClairS_TO_v0_3_0_pileup_nonsomatic", "DeepSomatic_TO_v1_8_0", "Longphase_TO_v0_0_1"]
    # run_calling_softwares = ["Longphase_TO_v0_0_1"]
    run_chr = None

    precision_calculator = PrecisionCalculator()
    commandrunner = CommandRunner()
    if("make" in get_config_column("make_command")):
        commandrunner.add_make_command()
        commandrunner.add_command(f"cp /bip8_disk/zhenyu112/longphase-to/longphase /bip8_disk/zhenyu112/longphase{tmp_num}")
    for variant_type in run_type:
        merge_list = []
        for sample in run_samples:
            sample_path = f"{data_root_path}/{sample}"
            sample_path_files = [d for d in os.listdir(sample_path)]
            # 選擇 ans
            ans_dir = max(sample_path_files, key=lambda file: ans_rank.get(file, -1)) if sample_path_files else "No valid files found"
            ans_path = f"{sample_path}/{ans_dir}"
            ans_vcf_name = [d for d in os.listdir(ans_path) if d.endswith('.vcf')]
            ans_vcf_name = [name for name in ans_vcf_name if all(filter not in name for filter in ans_vcf_filter)]
            ans_vcf_name = next((d for d in ans_vcf_name if variant_type in d.lower()), ans_vcf_name[0] if ans_vcf_name else None)
            ans_vcf_path = f"{ans_path}/{ans_vcf_name}"
            ans_bed = next((d for d in os.listdir(ans_path)if d.endswith('.bed')), None)
            ans_bed_path = f"{ans_path}/{ans_bed}"
            # 選擇 sequencing_platform
            sequencing_platforms = [file for file in sample_path_files if sequencing_platform_rank.get(file, -1) > 0]
            for sequencing_platform in sequencing_platforms:
                logging.info(f"{AnsiColors.GREEN}Runner:{sample} {sequencing_platform} {AnsiColors.RESET}")
                sequencing_platform_path = f"{sample_path}/{sequencing_platform}"
                subsample_path = f"{sequencing_platform_path}/subsample"
                if not os.path.exists(subsample_path):
                    logging.warning(f"{AnsiColors.RED}Subsample path does not exist: {subsample_path}. Skipping...{AnsiColors.RESET}")
                    continue
                sample_platform = f"{sample}_{sequencing_platform}"
                output_sample_platform_path = f"{output_root_path}/{sample_platform}"
                merge_list.append(sample_platform)
                purity_codes = [d for d in os.listdir(subsample_path)]
                puritys_dict = {p: v for p in purity_codes if (v := calculate_purity(p)) in run_puritys}
                if "ori" in run_puritys:
                    puritys_dict["ori"] = 1
                puritys_dict = dict(sorted(puritys_dict.items(), key=lambda item: item[1], reverse=True))

                for purity, value in puritys_dict.items():
                    name = purity
                    purity_path = ""
                    if name != "ori":
                        purity_path = f"{subsample_path}/{name}"
                    else:
                        purity_path = sequencing_platform_path
                    bam_file = next((f for f in os.listdir(purity_path) if f.endswith('.bam') and 'bl' not in f), None)
                    bam_path = f"{purity_path}/{bam_file}"
                    output_path = f"{output_sample_platform_path}/{name}/output"
                    for software in run_calling_softwares:
                        software_name, suffix = split_version_string(software)
                        calling_vcf_path = f"{purity_path}/{software_name}"
                        calling_vcf_df = None

                        match software_name:
                            case 'DeepSomatic_TO_v1_8_0':
                                calling_vcf_path += "/output.vcf"
                            case "ClairS_TO_v0_3_0":
                                if suffix == "":
                                    calling_vcf_path += f"/{variant_type}.vcf"
                                elif suffix == "_pileup":
                                    calling_vcf_path += f"/tmp/vcf_output/{variant_type}_pileup.vcf"
                                elif suffix == "_pileup_nonsomatic":
                                    calling_vcf_path += f"/tmp/vcf_output/{variant_type}_pileup_nonsomatic_tagging.vcf"
                            case "Longphase_TO_v0_0_1":
                                calling_vcf_path = f"{output_path}.vcf"
                                vcf_path = ""
                                if "DeepSomatic_TO_v1_8_0" in input_calling_softwares:
                                    deepvariant_vcf_path = f"{purity_path}/DeepSomatic_TO_v1_8_0/output.vcf"
                                    clairs_vcf_path = f"{purity_path}/ClairS_TO_v0_3_0/tmp/vcf_output/snv_pileup.vcf"
                                    vcf_path = f"{output_path}_merged.vcf"
                                    if not os.path.exists(output_path):
                                        os.makedirs(output_path)
                                    # subprocess.run(f"cp {deepvariant_vcf_path} {vcf_path}", shell=True, check=True)
                                    # subprocess.run(f'grep "^#"  {deepvariant_vcf_path} > {vcf_path}', shell=True, check=True)
                                    # subprocess.run(f'grep -v "^#" {deepvariant_vcf_path} | grep "PASS" | sed "s/1\\/1/0\\/1/g" >> {vcf_path}', shell=True, check=True)
                                    # subprocess.run(f'grep -v "^#" {deepvariant_vcf_path} | grep "GERMLINE" | sed "s/0\\/0/0\\/1/g" >> {vcf_path}', shell=True, check=True)
                                    # subprocess.run(f'grep -v "^#" {deepvariant_vcf_path} | grep "PON" | sed "s/1\\/1/0\\/1/g" >> {vcf_path}', shell=True, check=True)
                                    # subprocess.run(f'grep -v "^#" {clairs_vcf_path} | sed "s/:AF:/:VAF:/g" >> {vcf_path}', shell=True, check=True)
                                if "ClairS_TO_v0_3_0" in input_calling_softwares:
                                    vcf_path = f"{purity_path}/ClairS_TO_v0_3_0/tmp/vcf_output/snv_pileup.vcf"
                                    vcf_path = f"{purity_path}/ClairS_TO_v0_3_0/tmp/vcf_output/snv_pileup.vcf"
                                    if variant_type == "indel":
                                        snv_vcf_path = vcf_path
                                        indel_vcf_path = f"{purity_path}/ClairS_TO_v0_3_0/tmp/vcf_output/indel_pileup.vcf"
                                        vcf_path = f"{output_path}_merged.vcf"
                                        if not os.path.exists(output_path):
                                            os.makedirs(output_path)
                                        subprocess.run(f"cp {snv_vcf_path} {vcf_path}", shell=True, check=True)
                                        subprocess.run(f'grep -v "^#" {indel_vcf_path} >> {vcf_path}', shell=True, check=True)
                                    if not os.path.exists(vcf_path):
                                        logging.warning(f"{AnsiColors.RED}no found:{vcf_path}{AnsiColors.RESET}")
                                        continue
                                output_path +="2"
                                calling_vcf_path = f"{output_path}.vcf"
                                commandrunner.add_phase_command(bam_path, vcf_path, output_path, indel=variant_type == "indel", tmp_num=tmp_num)
                                commandrunner.add_merge_command(output_path, "somatic", run_mode="cat") # 合併
                                commandrunner.add_merge_command(output_path, "cluster", run_mode="cat") # 合併
                                # commandrunner.add_merge_command(output_path, "cnv", run_mode="cat") # 合併
                                commandrunner.add_merge_command(output_path, "purity", run_mode="cat") # 合併
                                commandrunner.add_vcf_index_command(output_path, background=True)
                                commandrunner.add_command(f"cp {output_path}.q1_q3 {output_sample_platform_path}/{name}.q1_q3")
                                # commandrunner.add_tag_command(bam_path, output_path, output_path)
                                # # commandrunner.add_tag_command(bam_path, output_path+"_tags", output_path)
                                # commandrunner.add_index_command(output_path+".bam")
                                # commandrunner.run(log_file_path=f"{output_root_path}/log/execution.log")
                                # # , test=True
                                # # , hide_log=True
                                # plot_vaf_distribution(f"{output_path}_cnv", f"{plot_vaf_distribution_path}/{sample_platform}/{purity_path}")

                                # tmp_path = f"{output_path}_somatic.vcf"
                                # subprocess.run(f"awk '/^#/ || /0\\|0/ || /\\.\\|0/ || /0\\|\\./' {calling_vcf_path} > {tmp_path}", shell=True, check=True)
                                # calling_vcf_path = tmp_path


                                # calling_vcf_df = BioinfoFileReader(calling_vcf_path).reader().cut_format_sample(["AF"]).df
                                test_path = f"/big8_disk/mingen112/test_data/HCC1954/ONT/subsample_bam/{name}/clairS/clairs-output-v0.4.0-ssrs/tmp/vcf_output/pileup.vcf"
                                # vcf_path = f"/big8_disk/data/HCC1954/ONT/subsample/t50_n00/ClairS_v0_4_0/output.vcf"
                                calling_vcf_df = BioinfoFileReader(test_path).reader().cut_format_sample(["AF"]).df.query("FILTER == 'PASS'")
                                calling_vcf_df["POS"] = calling_vcf_df["POS"] - 1
                                test = calling_vcf_df
                                tmp_name = "_ans_vaf"
                                # test = get_result(precision_calculator, f"{sample_platform},{value},no_cluster_qual",
                                #             in_vcf_path=calling_vcf_path, in_vcf_df=calling_vcf_df, ans_vcf_path=ans_vcf_path, ans_bed_path=ans_bed_path,
                                #             test_chr=run_chr, check_type=variant_type).query("ans == 'tp'")
                                
                                # # # # 创建AF值分布统计图表
                                if not test.empty and 'AF' in test.columns:
                                    # 确保输出目录存在
                                    os.makedirs(os.path.dirname(f"{output_sample_platform_path}/{name}_ans_vaf.png"), exist_ok=True)
                                    test['AF'] = test['AF'].astype(float)
                                    test['AF'] = test['AF'].round(2)
                                    # 统计AF值的分布
                                    af_counts = test['AF'].value_counts().sort_index()

                                    # 创建图表
                                    plt.figure(figsize=(12, 6)) # Increased figure size for better readability

                                    # 使用 plt.bar 繪製長條圖，設定寬度以避免重疊
                                    # Use a small width like 0.01 since AF values are close
                                    plt.bar(af_counts.index, af_counts.values, color='skyblue', width=0.008)

                                    plt.title('AF Distribution for True Positives') # More descriptive title
                                    plt.xlabel('Allele Frequency (AF)')
                                    plt.ylabel('Count')

                                    # 设置 x 轴范围和刻度
                                    plt.xlim(0, 1)
                                    x_ticks = np.arange(0, 1.1, 0.05) # Ticks from 0 to 1, step 0.05
                                    plt.xticks(x_ticks, [f"{tick:.2f}" for tick in x_ticks], rotation=45, ha='right') # Format labels to 2 decimal places

                                    plt.grid(True, axis='y', linestyle='--', alpha=0.7) # Add y-axis grid

                                    plt.tight_layout()

                                    # 保存图表
                                    plt.savefig(f"{output_sample_platform_path}/{name}_ans_vaf.png", dpi=300)
                                    plt.close()
                                    print(f"已生成AF值分布统计图表，保存在 {output_sample_platform_path}/{name}_ans_vaf.png")


                                # get_result(precision_calculator, f"{sample_platform},{value},total",
                                #             in_vcf_path=calling_vcf_path, in_vcf_df=None, ans_vcf_path=ans_vcf_path, ans_bed_path=ans_bed_path,
                                #             test_chr=run_chr, check_type=variant_type)

                                # calling_vcf_df = BioinfoFileReader(calling_vcf_path).reader().df
                                # calling_vcf_df["POS"] = calling_vcf_df["POS"] - 1
                                # calling_vcf_df["QUAL"] = calling_vcf_df["QUAL"].astype(float)
                                # cluster_df = BioinfoFileReader(output_path+"_cluster", columns=["CHROM", "POS", "count"]).reader().df
                                # calling_no_cluster_df = calling_vcf_df.merge(cluster_df, on=["CHROM", "POS"], how="left").query("count < 2 or count.isnull() or QUAL > 12")
                                # get_result(precision_calculator, f"{sample_platform},{value},no_cluster_qual",
                                #             in_vcf_path=calling_vcf_path, in_vcf_df=calling_no_cluster_df, ans_vcf_path=ans_vcf_path, ans_bed_path=ans_bed_path,
                                #             test_chr=run_chr, check_type=variant_type)
                                # calling_no_cluster_df = calling_vcf_df.merge(cluster_df, on=["CHROM", "POS"], how="left").query("count < 2 or count.isnull()")
                                # get_result(precision_calculator, f"{sample_platform},{value},no_cluster",
                                #             in_vcf_path=calling_vcf_path, in_vcf_df=calling_no_cluster_df, ans_vcf_path=ans_vcf_path, ans_bed_path=ans_bed_path,
                                #             test_chr=run_chr, check_type=variant_type)

                                # bed_df = BioinfoFileReader(output_path+"_SGE.bed", second_column=2).reader().df
                                # cluster_in_sge_df = Diff_bed_vcf(bed_df, cluster_df, 3).find_range_muti_chrs("SGE").vcf_df.query("SGE == 1")
                                # calling_no_cluster_in_sge_df = calling_vcf_df.merge(cluster_in_sge_df, on=["CHROM", "POS"], how="left").query("count < 2 or count.isnull()").drop(columns=["count"])
                                # calling_no_cluster_in_sge_df = calling_no_cluster_in_sge_df.merge(cluster_in_sge_df, on=["CHROM", "POS"], how="left").query("count < 2 or count.isnull() or QUAL > 12")
                                # get_result(precision_calculator, f"{sample_platform},{value},no_cluster_in_sge_qual",
                                #             in_vcf_path=calling_vcf_path, in_vcf_df=calling_no_cluster_in_sge_df, ans_vcf_path=ans_vcf_path, ans_bed_path=ans_bed_path,
                                #             test_chr=run_chr, check_type=variant_type)
                                # calling_no_cluster_in_sge_df = calling_vcf_df.merge(cluster_in_sge_df, on=["CHROM", "POS"], how="left").query("count < 2 or count.isnull()")
                                # get_result(precision_calculator, f"{sample_platform},{value},no_cluster_in_sge",
                                #             in_vcf_path=calling_vcf_path, in_vcf_df=calling_no_cluster_in_sge_df, ans_vcf_path=ans_vcf_path, ans_bed_path=ans_bed_path,
                                #             test_chr=run_chr, check_type=variant_type)

                                # commandrunner.add_vcf_index_command(calling_vcf_path.replace(".vcf", ""), background=True)
                                # commandrunner.run(log_file_path=f"{output_root_path}/log/execution.log")

                                # somatic_df = somatic_test(f"{output_path}_somatic")
                                # somatic_df['somatic'] = (somatic_df['HIGHV'] > 0) | ((somatic_df['LOWV'] / (somatic_df['LOWV'] + somatic_df['DISAGREEV'])) > 0.2)
                                # somatic_df.query("somatic == True", inplace=True)

                                # calling_ans_vcf_df = get_result(precision_calculator, f"{sample_platform},{value},{software},{variant_type}",
                                #                        in_vcf_path=calling_vcf_path, in_vcf_df=somatic_df, ans_vcf_path=ans_vcf_path, ans_bed_path=ans_bed_path,
                                #                        test_chr=run_chr, check_type=variant_type)
                                # calling_ans_vcf_df = calling_ans_vcf_df.query("somatic == True or ans != 'fp'")
                                # calling_ans_vcf_df = calling_ans_vcf_df.query("somatic == False and ans == 'tp'")
                                # calling_ans_vcf_df.loc[(calling_ans_vcf_df['somatic'] == False) & (calling_ans_vcf_df['ans'] == 'tp'), 'ans'] = 'fn'
                                # print(calling_ans_vcf_df)
                                # print(calling_ans_vcf_df['HIGHV'].value_counts())
                                # print(calling_ans_vcf_df['LOWV'].value_counts())
                                # calling_ans_vcf_df['test'] = calling_ans_vcf_df['HIGHV'] + calling_ans_vcf_df['LOWV']
                                # print(calling_ans_vcf_df['test'].value_counts())

                                # calling_ans_vcf_df['ratio'] = calling_ans_vcf_df['HIGHV'] / (calling_ans_vcf_df['HIGHV']+calling_ans_vcf_df['DISAGREEV'] + calling_ans_vcf_df['LOWV'])
                                # print(calling_ans_vcf_df['ratio'].value_counts())


                                # somatic_df = somatic_test(f"{output_path}_somatic").query("HIGHV == 0 and LOWV > 0")
                                # # somatic_df = somatic_test(f"{output_path}_somatic").query("HIGHV > 0 or ratio_t > 0.1")
                                # somatic_df = get_result(precision_calculator, f"{sample_platform},{value},{software}_diss3-add6,{variant_type}",
                                #             in_vcf_path=calling_vcf_path, in_vcf_df=somatic_df, ans_vcf_path=ans_vcf_path, ans_bed_path=ans_bed_path,
                                #             test_chr=run_chr, check_type=variant_type).query("ans == 'tp' or ans == 'fp'")
                                # tmp_name="_call"
                                # plot_ratio_t_distribution(somatic_df, f"{output_sample_platform_path}/{name}{tmp_name}.png")

                                # purity_df = purity_test(f"{output_path}_purity")
                                # purity_df = get_result(precision_calculator, f"{sample_platform},{value},{software},{variant_type}",
                                #                        in_vcf_path=calling_vcf_path, in_vcf_df=purity_df, ans_vcf_path=ans_vcf_path, ans_bed_path=ans_bed_path,
                                #                        test_chr=run_chr, check_type=variant_type).query("ans == 'tp'")
                                # purity_image_path = f"{output_sample_platform_path}/{name}.png"
                                # # purity_df = purity_df[purity_df['VAF'].notna() & (purity_df['VAF'] != 0)]
                                # purity_df = purity_df.merge(somatic_df, on=["CHROM", "POS"], how="left")
                                # # purity_df = purity_df[(purity_df["HIGHV"] ==0) &( purity_df["ratio_t"] > 0.2)]
                                # purity_df = purity_df[(purity_df["HIGHV"] > 0)]
                                # # purity_df = purity_df[(purity_df["HIGH"] == 0) | (purity_df["HIGH"].isna())]

                                # purity_df = purity_test(f"{output_path}_purity")
                                # purity_image_path = f"{output_sample_platform_path}/{name}.png"
                                # estimate_purity(purity_df, purity_image_path, name=f"{sample_platform},{value},{software},{variant_type}")
                                # calling_vcf_df = purity_df
                                # calling_ans_vcf_df = get_result(precision_calculator, f"{sample},{value},{software},{variant_type}",
                                #                     #    in_vcf_path=calling_vcf_path, in_vcf_df=calling_vcf_df, ans_vcf_path=ans_vcf_path, ans_bed_path=ans_bed_path,
                                #                        in_vcf_path=calling_vcf_path, in_vcf_df=calling_vcf_df, ans_vcf_path=ans_vcf_path, ans_bed_path=None,
                                #                        test_chr=run_chr, check_type=variant_type).query("ans == 'tp' or ans == 'fp'")
                                # estimate_purity(calling_ans_vcf_df, purity_image_path, name=f"{sample},{value},{software},{variant_type}")

                                # # calling_vcf_path = phased_tag_ans(output_path, ans_vcf_path, ans_bed_path)
                                # clairs_to_vcf_path = f"{purity_path}/ClairS_TO_v0_3_0/{variant_type}.vcf"
                                # reader = BioinfoFileReader(clairs_to_vcf_path).reader()
                                # reader.df = reader.df.query("FILTER == 'PASS'")
                                # reader.df["POS"] = reader.df["POS"] - 1
                                # calling_ans_vcf_df = get_result(precision_calculator, f"{sample},{value},ClairS_TO_v0_3_0,{variant_type}",
                                #                        in_vcf_path=clairs_to_vcf_path, in_vcf_df=reader.df, ans_vcf_path=ans_vcf_path, ans_bed_path=ans_bed_path,
                                #                        test_chr=run_chr, check_type=variant_type).query("ans == 'tp'")
                                # reader.df = calling_ans_vcf_df
                                # reader.cut_format_sample(['AF'])
                                # calling_ans_vcf_df = get_result(precision_calculator, f"{sample},{value},{software},{variant_type}",
                                #                        in_vcf_path=calling_vcf_path, in_vcf_df=None, ans_vcf_path=ans_vcf_path, ans_bed_path=ans_bed_path,
                                #                        test_chr=run_chr, check_type=variant_type).query("ans == 'tp'")
                                # unique_df = reader.df[~reader.df.set_index(['CHROM', 'POS']).index.isin(calling_ans_vcf_df.set_index(['CHROM', 'POS']).index)]
                                # # 轉成 float，並強制賦值回原欄位
                                # unique_df['AF'] = pd.to_numeric(unique_df['AF'], errors='coerce')
                                # # unique_df['AF_str'] = unique_df['AF'].apply(lambda x: format(float(x), '.2f'))
                                # unique_df.loc[:, 'AF'] = unique_df['AF'].round(2)  # 四捨五入到小數點第二位
                                # print(unique_df['AF'].value_counts().sort_index())

                                # # Create a bar plot for the value counts of unique_df['AF']
                                # af_counts = unique_df['AF'].value_counts().sort_index()
                                # plt.figure(figsize=(10, 6))
                                # af_counts.plot(kind='bar', color='skyblue')
                                # plt.title('Frequency of AF Values')
                                # plt.xlabel('AF Values')
                                # plt.ylabel('Frequency')
                                # plt.xticks(rotation=45)
                                # plt.tight_layout()
                                # plt.savefig(f"{output_path}_af.png")  # Save the figure
                                # plt.show()  # Display the figure


                                # purity 測試
                                # baf = np.array(purity_df['VAF'].tolist())
                                # test = create_distance_matrix(baf)
                                # plot_heatmap(test, f"{output_root_path}/{sample}/{name}_vote.png")
                                # print(test.shape)

                        if not calling_vcf_path or not os.path.exists(calling_vcf_path) and os.path.exists(f"{calling_vcf_path}.gz"):
                            decompress_gzip(f"{calling_vcf_path}.gz")
                        if not os.path.exists(calling_vcf_path):
                            logging.warning(f"{AnsiColors.RED}no found:{calling_vcf_path}{AnsiColors.RESET}")
                            # my_address, boos_address, boos_default_path = get_my_address(), get_boos_address(), get_boos_default_path()
                            # print(f"scp -q -J {my_address} -r {boos_address}:{boos_default_path}/{sample}/{purity}/{software} {data_root_path}/{sample}/{sequencing_platform}/subsample/{purity}/")
                            # subprocess.run(f"scp -q -J {my_address} -r {boos_address}:{boos_default_path}/{sample}/{purity}/{software} {data_root_path}/{sample}/{sequencing_platform}/subsample/{purity}/", shell=True, check=True)
                            continue
                        # test = get_result(precision_calculator, f"{sample_platform},{value},{software},{variant_type}",
                        #            in_vcf_path=calling_vcf_path, in_vcf_df=calling_vcf_df, ans_vcf_path=ans_vcf_path, ans_bed_path=ans_bed_path,
                        #            test_chr=run_chr, check_type=variant_type)
                        # print(test.query("ans == 'tp'"))
                        # print(test.query("ans == 'fp'"))
                        # get_result(precision_calculator, f"{sample_platform},{value},{software},{variant_type}",
                        #            in_vcf_path=calling_vcf_path, in_vcf_df=None, ans_vcf_path=ans_vcf_path, ans_bed_path=ans_bed_path,
                        #            test_chr=run_chr, check_type=variant_type)
                        # precision_calculator.print_confusion_matrix_metrics().write_results(f"{output_root_path}/metrics.txt")
                # merge_purity_images(f"{output_root_path}/{sample_platform}", f"{output_root_path}/{sample_platform}", list(puritys_dict.keys()), orientation='vertical', tmp_name = tmp_name)
                merge_purity_images(f"{output_root_path}/{sample_platform}", f"{output_root_path}/{sample_platform}", list(puritys_dict.keys()), orientation='horizontal', tmp_name = tmp_name)
                # puritys
                # logging.info(puritys_dict.keys())
            # sequencing_platforms
        # run_samples
        
        # merge_purity_images(f"{output_root_path}/merge", f"{output_root_path}", merge_list, orientation='horizontal', tmp_name = tmp_name)
        # merge_purity_images(f"{output_root_path}/merge", f"{output_root_path}", merge_list, orientation='vertical', tmp_name = tmp_name)
    # run_type
    # precision_calculator.print_confusion_matrix_metrics().write_results(f"{output_root_path}/metrics.txt")
    # precision_calculator.print_confusion_matrix_metrics().output_picture(f"{output_root_path}", format=["sample", "purity", "software", "variant_type"]).clear_results()

def calculate_loh_ratio():

    hg_tables_path = "/bip8_disk/zhenyu112/bioinfo_handler/data/hgTables.bed"
    hg_tables_df = BioinfoFileReader(hg_tables_path, second_column=1).reader().df
    hg_tables_df['POS1'] = hg_tables_df['POS1'].astype(int)  # 指定型態為 int32
    hg_tables_df['POS2'] = hg_tables_df['POS2'].astype(int)  # 指定型態為 int32
    hg_tables_df['length'] = hg_tables_df['POS2'] - hg_tables_df['POS1']
    total_length = hg_tables_df['length'].sum()
    print(total_length)

    dict_loh_ratio = {}
    one_loh_bed = f"{get_config_column('run_root_path')}/{get_config_column('output_folder')}"
    for sample in os.listdir(one_loh_bed):
        if sample == "log" or not os.path.isdir(f"{one_loh_bed}/{sample}"):
            continue  # Skip the log directory
        sample_path = f"{one_loh_bed}/{sample}"
        for purity in os.listdir(sample_path):
            purity_path = f"{sample_path}/{purity}"
            if purity == "log" or not os.path.isdir(f"{sample_path}/{purity}"):
                continue  # Skip the log directory
            bed_loh_df = BioinfoFileReader(purity_path+"/output_LOH.bed", second_column=2).reader().df
            bed_loh_df['POS1'] = bed_loh_df['POS1'].astype(int)
            bed_loh_df['POS2'] = bed_loh_df['POS2'].astype(int)
            bed_loh_df['length'] = bed_loh_df['POS2'] - bed_loh_df['POS1']
            loh_length = bed_loh_df['length'].sum()
            ratio = round(loh_length* 100/total_length, 2)
            dict_loh_ratio[f"{sample},{purity}"] = ratio
    for key, value in sorted(dict_loh_ratio.items(), key=lambda item: item[1], reverse=False):
        print(key,"\t", value)

def plot_ratio_t_distribution(somatic_df, output_path):
    """
    根据somatic_df中的ans列（tp和fp）分别生成两张图表，使用ratio_t作为x轴，统计ratio_t作为y轴
    
    参数:
        somatic_df: 包含ans列和ratio_t列的DataFrame
        output_path: 输出图片的路径
    """
    import matplotlib.pyplot as plt
    import numpy as np
    import os
    
    # 确保输出目录存在
    output_dir = os.path.dirname(output_path)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # 分离tp和fp数据
    tp_df = somatic_df[somatic_df['ans'] == 'tp']
    fp_df = somatic_df[somatic_df['ans'] == 'fp']
    
    # 创建两个子图
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    
    # 设置ratio_t的区间
    bins = np.linspace(0, 1, 21)  # 0到1之间创建20个区间
    
    # 绘制tp的直方图
    ax1.hist(tp_df['ratio_t'], bins=bins, alpha=0.7, color='blue', edgecolor='black')
    ax1.set_title('True Positives (TP)')
    ax1.set_xlabel('ratio_t')
    ax1.set_ylabel('Count')
    ax1.grid(True, linestyle='--', alpha=0.7)
    # ax1.set_ylim(0, 2500)  # 固定 y 轴高度
    
    # 绘制fp的直方图
    ax2.hist(fp_df['ratio_t'], bins=bins, alpha=0.7, color='red', edgecolor='black')
    ax2.set_title('False Positives (FP)')
    ax2.set_xlabel('ratio_t')
    ax2.set_ylabel('Count')
    ax2.grid(True, linestyle='--', alpha=0.7)
    # ax2.set_ylim(0, 2500)  # 固定 y 轴高度
    
    # 添加统计信息
    tp_mean = tp_df['ratio_t'].mean()
    tp_median = tp_df['ratio_t'].median()
    fp_mean = fp_df['ratio_t'].mean()
    fp_median = fp_df['ratio_t'].median()
    
    ax1.text(0.05, 0.95, f'Mean: {tp_mean:.3f}\nMedian: {tp_median:.3f}\nCount: {len(tp_df)}', 
             transform=ax1.transAxes, verticalalignment='top',
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    ax2.text(0.05, 0.95, f'Mean: {fp_mean:.3f}\nMedian: {fp_median:.3f}\nCount: {len(fp_df)}', 
             transform=ax2.transAxes, verticalalignment='top',
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    # 调整布局并保存图片
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()
    
    print(f"已生成ratio_t分布图，保存在 {output_path}")

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format='%(asctime)s  %(levelname)s: %(message)s')
    # logging.basicConfig(level=logging.DEBUG, format="%(asctime)s  %(levelname)s: %(message)s")
    file_for()
    # test_polynomial_regression()
    # calculate_loh_ratio()
    # test3()
    # test6()
    # merge_images()
    # slip_image()
    # gui_slip_image()
    # process_image_transform()
    # write_image()
    