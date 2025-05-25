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
    degree = 2
    # Get all subdirectories in root_dir
    all_dir = []
    check_dir = []
    
    # Scan all sample directories
    for sample in os.listdir(root_dir):
        sample_dir = os.path.join(root_dir, sample)
        if os.path.isdir(sample_dir) and not sample.endswith('log'):
        # if os.path.isdir(sample_dir) and not sample.endswith('log') and not sample.endswith('Dorado'):
            all_dir.append(sample_dir)
    print(all_dir)

    test_root_dir = f"{run_root_path}/{get_config_column('test_folder')}"
    # Scan all sample directories
    # for sample in os.listdir(test_root_dir):
    #     sample_dir = os.path.join(test_root_dir, sample)
    #     if os.path.isdir(sample_dir) and not sample.endswith('log'):
    #     # if os.path.isdir(sample_dir) and not sample.endswith('log') and not sample.endswith('Dorado'):
    #         # all_dir.append(sample_dir)
    #         check_dir.append(sample_dir)
    # print(check_dir)
    check_dir = all_dir
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
    # name_pattern = "*.csv"
    name_pattern = "*.q1_q3_tmp"
    name_pattern = "*.q1_q3"
    model_info = predictor.train(train_dir, model_path, name_pattern=name_pattern)
    # prediction = predictor.predict(check_dir, model_path, plot_model=True, use_function=False, name_pattern=name_pattern)
    # prediction = predictor.predict(check_dir, model_path, plot_model=True, use_function=True, name_pattern=name_pattern)
    # prediction = predictor.predict(check_dir, model_path, plot_model=False, use_function=True, name_pattern=name_pattern)
    prediction = predictor.predict(check_dir, model_path, plot_model=False, use_function=False, name_pattern=name_pattern)
    print_grouped_data(prediction)
    if isinstance(prediction, list) and len(prediction) > 0 and isinstance(prediction[0], tuple):
        plot_prediction_vs_actual(prediction, output_path=f"{root_dir}/purity_prediction_results.png")
        
    test_purity_table = """
purity	1	0.8	0.6	0.4	0.2
COLO829_ONT_PAO	0.99	0.76	0.54	0.4	0.42
H1437_ONT	1	0.79	0.6	0.41	1
H2009_ONT	0.99	0.83	0.65	0.45	0.32
HCC1395_ONT	0.99	0.73	0.6	0.49	1
HCC1395_ONT_Dorado	0.98	0.79	0.65	0.49	1
HCC1937_ONT	0.99	0.78	0.59	1	1
HCC1954_ONT	null	null	null	1	1
"""
    test_purity_table = """
purity	1	0.8	0.6	0.4	0.2
COLO829_ONT_PAO	0.74	0.63	0.42	0.23	null
H1437_ONT	0.84	0.69	0.47	0.29	0.24
H2009_ONT	0.92	0.73	0.51	0.32	0.25
HCC1395_ONT	0.9	0.68	0.46	0.4	null
HCC1395_ONT_Dorado	0.91	0.7	0.47	0.4	0.85
HCC1937_ONT	0.81	0.71	0.47	0.39	0.26
HCC1954_ONT	0.5	0.56	0.32	0.38	0.21
"""
    # purity_df = pd.read_csv(StringIO(test_purity_table), sep="\t")
    # # 取得純度欄位名稱（去掉第一欄 'purity'）
    # purity_levels = [float(x) for x in purity_df.columns[1:]]
    # # 組成 prediction-like 結構
    # prediction_from_table = []
    # for idx, row in purity_df.iterrows():
    #     sample_name = row['purity']
    #     for i, purity in enumerate(purity_levels):
    #         val = row.iloc[i+1]
    #         # 跳過 null 或 nan
    #         if isinstance(val, str) and val.lower() == 'null':
    #             continue
    #         # if isinstance(val, float) and (math.isnan(val) or val is None):
    #         #     continue
    #         actual_purity = purity
    #         predicted_purity = float(val)
    #         error = abs(predicted_purity - actual_purity)
    #         prediction_from_table.append((actual_purity, predicted_purity, error, sample_name))
    # print(prediction_from_table)
    # plot_prediction_vs_actual(prediction_from_table, output_path=f"{root_dir}/purity_prediction_results2.png")

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

def get_ref_length():
    hg_tables_path = "/bip8_disk/zhenyu112/bioinfo_handler/data/hgTables.bed"
    hg_tables_df = BioinfoFileReader(hg_tables_path, second_column=1).reader().df
    hg_tables_df['POS1'] = hg_tables_df['POS1'].astype(int)  # 指定型態為 int32
    hg_tables_df['POS2'] = hg_tables_df['POS2'].astype(int)  # 指定型態為 int32
    hg_tables_df['length'] = hg_tables_df['POS2'] - hg_tables_df['POS1']
    total_length = hg_tables_df['length'].sum()
    return total_length

def file_for():
    # ref_length = get_ref_length()
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
    plot_low_vote_distribution_path = f"{get_config_column('run_root_path')}/{get_config_column('output_plot_folder')}/low_vote"
    tmp_num = get_config_column("tmp")
    # run_type = ["snv"]  # snv indel
    # run_samples = ["H2009", "HCC1954", "HCC1937", "COLO829", "H1437", "HCC1395"]
    # run_samples = ["HCC1395"]
    # run_puritys = [1, 0.8, 0.6, 0.4, 0.2]
    # run_calling_softwares = ["ClairS_TO_v0_3_0", "ClairS_TO_v0_3_0_pileup", "ClairS_TO_v0_3_0_pileup_nonsomatic", "DeepSomatic_TO_v1_8_0", "Longphase_TO_ssrs_v0_0_1"]
    # run_calling_softwares = ["Longphase_TO_ssrs_v0_0_1"]
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
                    for software in run_calling_softwares:
                        output_path = f"{output_sample_platform_path}/{name}/output"
                        software_name, suffix = split_version_string(software)
                        calling_vcf_path = f"{purity_path}/{software_name}"
                        calling_vcf_df = None
                        if not os.path.exists(output_path):
                            output_dir = os.path.dirname(output_path)
                            os.makedirs(output_dir, exist_ok=True)
                        match software_name:
                            case 'DeepSomatic_TO_v1_8_0':
                                calling_vcf_path += "/output.vcf"
                            case "ClairS_TO_v0_3_0" | "ClairS_TO_ss_v0_3_0":
                                if suffix == "":
                                    calling_vcf_path += f"/{variant_type}.vcf"
                                elif suffix == "_pileup":
                                    calling_vcf_path += f"/tmp/vcf_output/{variant_type}_pileup.vcf"
                                elif suffix == "_pileup_nonsomatic":
                                    calling_vcf_path += f"/tmp/vcf_output/{variant_type}_pileup_nonsomatic_tagging.vcf"
                            case "Longphase_TO_ssrs_v0_0_1" | "Longphase_TO_ss_v0_0_1" | "Longphase_TO_deep_v0_0_1":
                                # calling_vcf_path = f"{output_path}.vcf"
                                vcf_path = ""
                                caller = ""
                                if "Longphase_TO_deep_v0_0_1" == software_name:
                                    output_path = f"{output_path}_deepsomatic_to"
                                    caller = "deepsomatic_to"
                                    vcf_path = f"{purity_path}/DeepSomatic_TO_v1_8_0/output.vcf"
                                if "Longphase_TO_ssrs_v0_0_1" == software_name or "Longphase_TO_ss_v0_0_1" == software_name:
                                    if "Longphase_TO_ss_v0_0_1" == software_name:
                                        output_path = f"{output_path}_clairs_to_ss"
                                        caller = "clairs_to_ss"
                                        input_calling_softwares = "ClairS_TO_ss_v0_3_0"
                                    elif "Longphase_TO_ssrs_v0_0_1" == software_name:
                                        output_path = f"{output_path}_clairs_to_ssrs"
                                        caller = "clairs_to_ssrs"
                                        input_calling_softwares = "ClairS_TO_v0_3_0"
                                    if "pileup" in input_calling_softwares:
                                        input_calling_softwares = input_calling_softwares.replace("_pileup", "")
                                        vcf_path = f"{purity_path}/{input_calling_softwares}/tmp/vcf_output/snv_pileup.vcf"
                                    else:
                                        vcf_path = f"{purity_path}/{input_calling_softwares}/snv.vcf"
                                    # vcf_path = f"/big8_disk/zhenyu112/clair_purity/{sample_lower}/{value}/ClairS_TO_v0_3_0ssrs/tmp/vcf_output/{variant_type}_pileup.vcf"
                                    if variant_type == "indel":
                                        snv_vcf_path = vcf_path
                                        if "pileup" in input_calling_softwares:
                                            indel_vcf_path = f"{purity_path}/{input_calling_softwares}/tmp/vcf_output/indel_pileup.vcf"
                                        else:
                                            indel_vcf_path = f"{purity_path}/{input_calling_softwares}/indel.vcf"
                                        vcf_path = f"{output_path}_merged.vcf"
                                        if not os.path.exists(output_path):
                                            os.makedirs(output_path)
                                        subprocess.run(f"cp {snv_vcf_path} {vcf_path}", shell=True, check=True)
                                        subprocess.run(f'grep -v "^#" {indel_vcf_path} >> {vcf_path}', shell=True, check=True)
                                    if not os.path.exists(vcf_path):
                                        logging.warning(f"{AnsiColors.RED}no found:{vcf_path}{AnsiColors.RESET}")
                                        continue

                                calling_vcf_path = f"{output_path}.vcf"
                                commandrunner.add_phase_command(bam_path, vcf_path, output_path, indel=variant_type == "indel", caller=caller, tmp_num=tmp_num)
                                # commandrunner.add_merge_command(output_path, "somatic", run_mode="cat") # 合併
                                # commandrunner.add_merge_command(output_path, "vote", run_mode="cat") # 合併
                                # commandrunner.add_merge_command(output_path, "cluster", run_mode="cat") # 合併
                                # commandrunner.add_merge_command(output_path, "cnv", run_mode="cat") # 合併
                                # commandrunner.add_merge_command(output_path, "purity", run_mode="cat") # 合併
                                # commandrunner.add_vcf_index_command(output_path, background=True)
                                commandrunner.add_command(f"cp {output_path}.q1_q3 {output_sample_platform_path}/{name}.q1_q3")
                                # commandrunner.add_tag_command(bam_path, output_path, output_path)
                                # # commandrunner.add_tag_command(bam_path, output_path+"_tags", output_path)
                                # commandrunner.add_index_command(output_path+".bam")
                                # commandrunner.run(log_file_path=f"{output_root_path}/log/execution.log")
                                # , test=True
                                # , hide_log=True
                                # # plot_vaf_distribution(f"{output_path}_cnv", f"{plot_vaf_distribution_path}/{sample_platform}/{purity_path}")

                                tmp_path = f"{output_path}_somatic.vcf"
                                # subprocess.run(f"awk '/^#/ || /0\\|0/ || /\\.\\|0/ || /0\\|\\./' {calling_vcf_path} > {tmp_path}", shell=True, check=True)
                                calling_vcf_path = tmp_path

                                # 排除小sge
                                # sge_df = BioinfoFileReader(f"{output_path}_SGE.bed",second_column= 2).reader().df
                                # somatic_df = BioinfoFileReader(calling_vcf_path).reader().df
                                # out_seg_somatic_df = Diff_bed_vcf(sge_df, somatic_df,1).find_range_muti_chrs("SGE").vcf_df
                                # out_seg_somatic_df["POS"] = out_seg_somatic_df["POS"] - 1
                                # ans_somatic_df = get_result(precision_calculator, f"{sample_platform},{value},{software},{variant_type}",
                                #            in_vcf_path=calling_vcf_path, in_vcf_df=out_seg_somatic_df, ans_vcf_path=ans_vcf_path, ans_bed_path=ans_bed_path,
                                #            test_chr=run_chr, check_type=variant_type).query("ans == 'tp' or ans == 'fp'")
                                # print(ans_somatic_df.value_counts(["ans","SGE"]))
                                # somatic_df["POS"] = somatic_df["POS"] - 1
                                # ans_somatic_df = get_result(precision_calculator, f"{sample_platform},{value},{software}_test,{variant_type}",
                                #            in_vcf_path=calling_vcf_path, in_vcf_df=out_seg_somatic_df.query("SGE == 0"), ans_vcf_path=ans_vcf_path, ans_bed_path=ans_bed_path,
                                #            test_chr=run_chr, check_type=variant_type).query("ans == 'tp' or ans == 'fp'")
                                # precision_calculator.print_confusion_matrix_metrics()

                                # 读取投票结果文件
                                # # vote_columns = ["CHROM", "POS", 
                                # #     "MID_HIGH_SR", "MID_HIGH_SA", "MID_HIGH_SR_VOTE", "MID_HIGH_SA_VOTE", 
                                # #     "LEFT_HIGH_SR", "LEFT_HIGH_SA", "LEFT_HIGH_SR_VOTE", "LEFT_HIGH_SA_VOTE", 
                                # #     "RIGHT_HIGH_SR", "RIGHT_HIGH_SA", "RIGHT_HIGH_SR_VOTE", "RIGHT_HIGH_SA_VOTE", 
                                # #     "MID_LOW", "LEFT_LOW", "RIGHT_LOW", "DISAGREE","VOTE", "VOTE_RAITO"]
                                # # vote_df = BioinfoFileReader(f"{output_path}_vote", columns=vote_columns).reader().df
                                # vote_columns = ["CHROM", "POS","HIGH", "LOW", "DISAGREE","TOTAL_RATIO"]
                                # vote_df = BioinfoFileReader(f"{output_path}_somatic", columns=vote_columns).reader().df
                                # vote_df["VOTE_RAITO"] = vote_df["TOTAL_RATIO"]/vote_df["HIGH"]
                                
                                # # # # 计算高置信度和低置信度的总和
                                # # high_cols = ["MID_HIGH_SR", "MID_HIGH_SA", "LEFT_HIGH_SR", "LEFT_HIGH_SA", 
                                # #             "RIGHT_HIGH_SR", "RIGHT_HIGH_SA", "MID_HIGH_SR_VOTE", "MID_HIGH_SA_VOTE",
                                # #             "LEFT_HIGH_SR_VOTE", "LEFT_HIGH_SA_VOTE", "RIGHT_HIGH_SR_VOTE", "RIGHT_HIGH_SA_VOTE"]
                                # # low_cols = ["MID_LOW", "LEFT_LOW", "RIGHT_LOW"]
                                
                                # # vote_df["total_high"] = vote_df[high_cols].sum(axis=1)
                                # # vote_df["total_low"] = vote_df[low_cols].sum(axis=1)
                                # # vote_df["low_ratio"] = vote_df["total_low"] / (vote_df["total_low"]+vote_df["DISAGREE"])
                                # vote_df["total_high"] = vote_df["HIGH"]
                                # vote_df["total_low"] = vote_df["LOW"]
                                # vote_df["low_ratio"] = vote_df["LOW"] / (vote_df["LOW"]+vote_df["DISAGREE"])
                                # print(vote_df.shape)

                                # vote_df.query("total_high >= 1 or low_ratio >= 0.2", inplace=True)
                                # print(vote_df.shape)
                                # # vote_df.query("total_high < 1 and low_ratio >= 0.2", inplace=True)
                                # # vote_df.query("total_high >= 1 ", inplace=True)
                                # ori_vcf_df = get_result(precision_calculator, f"{sample_platform},{value},{software}_ss,{variant_type}",
                                #             in_vcf_path=calling_vcf_path, in_vcf_df=vote_df, ans_vcf_path=ans_vcf_path, ans_bed_path=ans_bed_path,
                                #             test_chr=run_chr, check_type=variant_type).query("ans == 'tp' or ans == 'fp'")
                                # vote_df["POS"] = vote_df["POS"] - 1
                                # # print(ori_vcf_df.query("ans == 'tp' and total_high > 20 and  VOTE_RAITO < 0.8"))
                                # # print(ori_vcf_df.query("ans == 'tp' and  VOTE_RAITO < 0.8"))
                                # # print(ori_vcf_df.query("ans == 'fp' and  VOTE_RAITO < 0.8"))
                                # print(ori_vcf_df.query("ans == 'fp'"))

                                # plot_ratio_t_distribution(ori_vcf_df, f"{output_sample_platform_path}/{name}_ss_votes_ratio.png", column_name="VOTE_RAITO")
                                # tmp_name = "_ss_votes_ratio"

                                # vote_df.query("(total_high >= 1 and VOTE_RAITO > 0.8) or low_ratio >= 0.2", inplace=True)
                                # ori_vcf_df = get_result(precision_calculator, f"{sample_platform},{value},{software}_ss_new,{variant_type}",
                                #             in_vcf_path=calling_vcf_path, in_vcf_df=vote_df, ans_vcf_path=ans_vcf_path, ans_bed_path=ans_bed_path,
                                #             test_chr=run_chr, check_type=variant_type).query("ans == 'tp' or ans == 'fp'")
                                # precision_calculator.print_confusion_matrix_metrics()

                                # test_df = get_result(precision_calculator, f"{sample_platform},{value},{software}test,{variant_type}",
                                #             in_vcf_path=calling_vcf_path, in_vcf_df=None, ans_vcf_path=ans_vcf_path, ans_bed_path=ans_bed_path,
                                #             test_chr=run_chr, check_type=variant_type).query("ans == 'tp' or ans == 'fp'")
                                # ori_vcf_df = ori_vcf_df.merge(test_df[["CHROM", "POS"]], on=["CHROM", "POS"], how="left", indicator=True)
                                # print(ori_vcf_df.drop(high_cols, axis=1).query("ans == 'tp' and _merge == 'left_only'"))

                                # plot_ssrs_votes_comparison(vote_df, high_cols, low_cols, precision_calculator,
                                #                            f"{sample_platform},{value},{software}_ssrs,{variant_type}", 
                                #                            variant_type, calling_vcf_path, ans_vcf_path, ans_bed_path, run_chr, 
                                #                            f"{output_sample_platform_path}/{name}_ssrs_votes_bar_one.png")

                                # output_path = output_path.replace("test-somatic-ss-log", "test-somatic-log")
                                # output_sample_platform_path = output_sample_platform_path.replace("test-somatic-ss-log", "test-somatic-log")
                                # vote_df = BioinfoFileReader(f"{output_path}_vote", columns=vote_columns).reader().df

                                
                                # plot_ssrs_votes_comparison(vote_df, high_cols, low_cols, precision_calculator,
                                #                            f"{sample_platform},{value},{software}_ss,{variant_type}", 
                                #                            variant_type, calling_vcf_path, ans_vcf_path, ans_bed_path, run_chr, 
                                #                            f"{output_sample_platform_path}/{name}_ss_votes_bar_one.png")
                                

                                # peileup_df = BioinfoFileReader(vcf_path).reader().df
                                # somatic_df = somatic_test(f"{output_path}_somatic")
                                # somatic_df['H'] = (somatic_df['HIGHV'] > 0)
                                # somatic_df['Lratio'] = (somatic_df['LOWV'] / (somatic_df['LOWV'] + somatic_df['DISAGREEV'])).round(2)
                                # somatic_df['L'] = ((somatic_df['LOWV'] > 0) & (somatic_df['H'] == False))
                                # # somatic_df['L'] = ((somatic_df['Lratio'] > 0.2) & (somatic_df['H'] == False))

                                # somatic_df.query("H == True or L == True", inplace=True)
                                # ori_vcf_df = get_result(precision_calculator, f"{sample_platform},{value},{software}_ss,{variant_type}",
                                #             in_vcf_path=calling_vcf_path, in_vcf_df=somatic_df, ans_vcf_path=ans_vcf_path, ans_bed_path=ans_bed_path,
                                #             test_chr=run_chr, check_type=variant_type).query("ans == 'tp' or ans == 'fp'")
                                # # precision_calculator.print_confusion_matrix_metrics().clear_results()
                                # print("ss calling total",ori_vcf_df.shape)

                                # output_high_path_picture = f"{output_sample_platform_path}/{name}_highv_distribution.png"
                                # output_low_path_picture = f"{output_sample_platform_path}/{name}_lowv_distribution.png"
                                # # plot_ratio_t_distribution(ori_vcf_df[ori_vcf_df['H'] == True], output_high_path_picture, column_name="HIGHV")
                                # # plot_ratio_t_distribution(ori_vcf_df[ori_vcf_df['L'] == True], output_low_path_picture, column_name="Lratio")

                                # old_output_path = output_path.replace("test-somatic-ss", "test-somatic")
                                # somatic_old_df = somatic_test(f"{old_output_path}_somatic")
                                # somatic_old_df.columns = somatic_old_df.columns.str.replace('HIGHV', 'HIGHVO')
                                # somatic_old_df.columns = somatic_old_df.columns.str.replace('LOWV', 'LOWVO')
                                # somatic_old_df.columns = somatic_old_df.columns.str.replace('DISAGREEV', 'DISAGREEVO')
                                # somatic_old_df['HO'] = (somatic_old_df['HIGHVO'] > 0)
                                # somatic_old_df['LratioO'] = (somatic_old_df['LOWVO'] / (somatic_old_df['LOWVO'] + somatic_old_df['DISAGREEVO'])).round(2)
                                # somatic_old_df['LO'] = ((somatic_old_df['LOWVO'] > 0) & (somatic_old_df['HO'] == False))

                                # somatic_old_df.query("HO == True or LO == True", inplace=True)
                                # ori_old_vcf_df = get_result(precision_calculator, f"{sample_platform},{value},{software}_ssrs,{variant_type}",
                                #             in_vcf_path=calling_vcf_path, in_vcf_df=somatic_old_df, ans_vcf_path=ans_vcf_path, ans_bed_path=ans_bed_path,
                                #             test_chr=run_chr, check_type=variant_type).query("ans == 'tp' or ans == 'fp'")
                                # print("ssrs calling total",ori_old_vcf_df.shape)
                                # # precision_calculator.print_confusion_matrix_metrics().clear_results()
                                # output_high_path_picture = f"{output_sample_platform_path}/{name}_highvo_distribution.png"
                                # output_low_path_picture = f"{output_sample_platform_path}/{name}_lowvo_distribution.png"
                                # # plot_ratio_t_distribution(ori_old_vcf_df[ori_old_vcf_df['HO'] == True], output_high_path_picture, column_name="HIGHVO")
                                # # plot_ratio_t_distribution(ori_old_vcf_df[ori_old_vcf_df['LO'] == True], output_low_path_picture, column_name="LratioO")

                                # # new_fp = ori_vcf_df.query("ans == 'fp' and HIGHV >= 1")
                                # # old_fp = ori_old_vcf_df.query("ans == 'fp' and HIGHVO >= 1")
                                # new_fp = ori_vcf_df.query("ans == 'fp'")
                                # old_fp = ori_old_vcf_df.query("ans == 'fp'")
                                # merged = new_fp.merge(
                                #     old_fp[["CHROM", "POS"]],
                                #     on=["CHROM", "POS"],
                                #     how="left",
                                #     indicator=True
                                # )
                                # existing_count = merged[merged["_merge"] == "both"].shape[0]
                                # non_existing_count = merged[merged["_merge"] == "left_only"].shape[0]
                                # non_existing_both_count = merged[merged["_merge"] == "both"].merge(somatic_old_df, on=["CHROM", "POS"], how="inner").shape[0]
                                # print(f"存在的數量: {existing_count}, 不存在的數量: {non_existing_count}, 不存在的數量且在old_vcf_df中的數量: {non_existing_both_count}")
                                # only_in_new = merged[merged["_merge"] == "left_only"].drop(columns=["_merge"])
                                # only_in_new.to_csv(f"only_in_new_fp.csv", index=False, sep=' ')

                                # test = peileup_df.merge(only_in_new, on=["CHROM", "POS"], how="inner")[merged["_merge"] == "both"]
                                # print(test)

                                # new_fp = ori_vcf_df.query("ans == 'tp'")
                                # old_fp = ori_old_vcf_df.query("ans == 'tp'")
                                # # # new_fp = ori_vcf_df.query("ans == 'tp' and HIGHV >= 1")
                                # # # old_fp = ori_old_vcf_df.query("ans == 'tp' and HIGHVO >= 1")
                                # merged = new_fp.merge(
                                #     old_fp[["CHROM", "POS"]],
                                #     on=["CHROM", "POS"],
                                #     how="left",
                                #     indicator=True
                                # )
                                # existing_count = merged[merged["_merge"] == "both"].shape[0]
                                # non_existing_count = merged[merged["_merge"] == "left_only"].shape[0]
                                # non_existing_both_count = merged[merged["_merge"] == "both"].merge(somatic_old_df, on=["CHROM", "POS"], how="inner").shape[0]
                                # print(f"存在的數量: {existing_count}, 不存在的數量: {non_existing_count}, 不存在的數量且在old_vcf_df中的數量: {non_existing_both_count}")
                                # only_in_new = merged[merged["_merge"] == "left_only"].drop(columns=["_merge"])
                                # only_in_new.to_csv(f"only_in_new_tp.csv", index=False, sep=' ')
                                
                                # test = peileup_df.merge(only_in_new, on=["CHROM", "POS"], how="inner")[merged["_merge"] == "both"]
                                # print(test)



                                # 讀取AF
                                # calling_vcf_df = BioinfoFileReader(calling_vcf_path).reader().cut_format_sample(["AF"]).df
                                # calling_vcf_df["ratio"] = calling_vcf_df["AF"].astype(float)

                                # ============================輸出purity圖片
                                # purity_df = purity_test(f"{output_path}_purity")
                                # purity_image_path = f"{output_sample_platform_path}/{name}.png"
                                # tmp_name=""
                                # estimate_purity(purity_df, purity_image_path, name=f"{sample_platform},{value},{software},{variant_type}")

                                # 計算LOH比例
                                # loh_df = BioinfoFileReader(f"{output_path}_LOH.bed").reader().df
                                # loh_df["length"] = loh_df["POS2"] - loh_df["POS1"]
                                # loh_length = loh_df["length"].sum()
                                # loh_ratio = loh_length / ref_length
                                # # loh_bool = 1 if loh_ratio > 0.1 else 0
                                # subprocess.run(f"echo 'lohratio,{loh_ratio:.4f}' >> {output_sample_platform_path}/{name}_1.csv", shell=True, check=True)
                                # with open(f"{output_sample_platform_path}/{name}.q1_q3_tmp", "w") as tmp_file:
                                #     first_line = open(f"{output_sample_platform_path}/{name}.q1_q3").readline()
                                #     tmp_file.write(first_line.strip())
                                #     tmp_file.write(f'\t{loh_ratio:.4f}')
                                #     tmp_file.write("\n")

                                # # # # 创建AF值分布统计图表
                                # if not test.empty and 'AF' in test.columns:
                                #     plot_af_distribution(test, f"{output_sample_platform_path}/{name}", file_suffix=f"_ans_vaf")

                                # purity 測試
                                # baf = np.array(purity_df['VAF'].tolist())
                                # test = create_distance_matrix(baf)
                                # plot_heatmap(test, f"{output_root_path}/{sample}/{name}_vote.png")
                                # print(test.shape)
                            case "Ascat":
                                normal_dir = f"{data_root_path}/{sample}/{sequencing_platform}/subsample/t00_n25/"
                                bam_files = [f for f in os.listdir(normal_dir) if f.endswith('.bam')][0]
                                normal_bam_path = normal_dir+bam_files
                                output_dir = f"{output_sample_platform_path}/{name}"
                                if sample == "COLO829" or sample == "H1437" or sample == "H2009":
                                    gender = "XY"
                                else:
                                    gender = "XX"
                                subprocess.run(f"bash /bip8_disk/zhenyu112/bioinfo_handler/scripts/run_ascat.sh {bam_path} {normal_bam_path} {output_dir} {gender}", shell=True, check=True)
                            case "ClairS_TO_Ascat":
                                command = f"cp -r {purity_path}/{input_calling_softwares}/* {output_sample_platform_path}/{name}/"
                                print(command)
                                subprocess.run(command, shell=True, check=True)
                                with open(f"{output_path}_log.txt", "w") as log_file:
                                    command = f"/bip8_disk/zhenyu112/ClairS-TO/run_clairs_to --tumor_bam_fn {bam_path} --ref_fn /mnt/ramdisk/GRCh38_no_alt_analysis_set.fasta --threads 64 --platform ont_r10_dorado_sup_5khz_ssrs --output_dir {output_sample_platform_path}/{name}"
                                    print(f"Running: {command}")
                                    process = subprocess.Popen(
                                        command,
                                        shell=True,
                                        stdout=log_file,
                                        stderr=log_file
                                    )
                                    process.communicate()
                                print(name)
                                subprocess.run(f"grep ================================================================= {output_path}_log.txt", shell=True, check=False)

                        if not calling_vcf_path or not os.path.exists(calling_vcf_path) and os.path.exists(f"{calling_vcf_path}.gz"):
                            decompress_gzip(f"{calling_vcf_path}.gz")
                        # subprocess.run(f"rm -r {data_root_path}/{sample}/{sequencing_platform}/subsample/{purity}/ClairS_TO_ss_v0_3_0", shell=True, check=True)

                        # print(f"ln-s /big8_disk/mingen112/test_data/COLO829_R10/ONT/subsample_bam/{purity}/COLO829_R10_{purity}.bam {data_root_path}/COLO829/ONT_R10/subsample/{purity}/")
                        # print(f"ln-s /big8_disk/mingen112/test_data/COLO829_R10/ONT/subsample_bam/{purity}/COLO829_R10_{purity}.bam.bai {data_root_path}/COLO829/ONT_R10/subsample/{purity}/")
                        # subprocess.run(f"ls -lh /big8_disk/mingen112/test_data/COLO829_R10/ONT/subsample_bam/{purity}", shell=True, check=True)
                        # # subprocess.run(f"mkdir -p {data_root_path}/COLO829/ONT_R10/subsample/{purity}/", shell=True, check=True)
                        # subprocess.run(f"ln -s /big8_disk/mingen112/test_data/COLO829_R10/ONT/subsample_bam/{purity}/COLO829_R10_{purity}.bam {data_root_path}/COLO829/ONT_R10/subsample/{purity}/", shell=True, check=True)
                        # subprocess.run(f"ln -s /big8_disk/mingen112/test_data/COLO829_R10/ONT/subsample_bam/{purity}/COLO829_R10_{purity}.bam.bai {data_root_path}/COLO829/ONT_R10/subsample/{purity}/", shell=True, check=True)
                        # subprocess.run(f"mv {data_root_path}/COLO829/ONT_R10/subsample/{purity}/COLO829_R10_{purity}.bam.bai {data_root_path}/COLO829/ONT_R10/subsample/{purity}/COLO829_{purity}.bam.bai", shell=True, check=True)
                        # subprocess.run(f"mv {data_root_path}/COLO829/ONT_R10/subsample/{purity}/COLO829_R10_{purity}.bam {data_root_path}/COLO829/ONT_R10/subsample/{purity}/COLO829_{purity}.bam", shell=True, check=True)
                        # subprocess.run(f"ls {data_root_path}/COLO829/ONT_R10/subsample/{purity}/", shell=True, check=True)
                        # subprocess.run(f"mv {data_root_path}/{sample}/{sequencing_platform}/subsample/t00_n25/{sample}_normal_25x.bam {data_root_path}/{sample}/{sequencing_platform}/subsample/t00_n25/{sample}_t00_n25.bam ", shell=True, check=True)
                        # subprocess.run(f"mv {data_root_path}/{sample}/{sequencing_platform}/subsample/t00_n25/{sample}_normal_25x.bam.bai {data_root_path}/{sample}/{sequencing_platform}/subsample/t00_n25/{sample}_t00_n25.bam.bai ", shell=True, check=True)
                        # subprocess.run(f"ls {data_root_path}/{sample}/{sequencing_platform}/subsample/t00_n25/", shell=True, check=True)

                        if not os.path.exists(calling_vcf_path):
                            logging.warning(f"{AnsiColors.RED}no found:{calling_vcf_path}{AnsiColors.RESET}")
                            my_address, boos_address, boos_default_path = get_my_address(), get_boos_address(), get_boos_default_path()
                            print(f"scp -q -J {my_address} -r {boos_address}:{boos_default_path}/{sample}/{sequencing_platform}/subsample/{purity}/{software} {data_root_path}/{sample}/{sequencing_platform}/subsample/{purity}/")
                            # subprocess.run(f"scp -q -J {my_address} -r {boos_address}:{boos_default_path}/{sample}/{sequencing_platform}/subsample/{purity}/{software} {data_root_path}/{sample}/{sequencing_platform}/subsample/{purity}/", shell=True, check=True)
                            continue
                        # get_result(precision_calculator, f"{sample_platform},{value},{software},{variant_type}",
                        #            in_vcf_path=calling_vcf_path, in_vcf_df=None, ans_vcf_path=ans_vcf_path, ans_bed_path=ans_bed_path,
                        #            test_chr=run_chr, check_type=variant_type)
                        no_filter = True
                        if "Deep" in software or "Clair" in software:
                            no_filter = False
                        get_result(precision_calculator, f"{sample_platform},{value},{software},{variant_type}",
                                   in_vcf_path=calling_vcf_path, in_vcf_df=None, ans_vcf_path=ans_vcf_path, ans_bed_path=ans_bed_path,
                                   test_chr=run_chr, check_type=variant_type, no_filter=no_filter)
                        precision_calculator.print_confusion_matrix_metrics().write_results(f"{output_root_path}/metrics_l.txt")
                # merge_purity_images(f"{output_root_path}/{sample_platform}", f"{output_root_path}/{sample_platform}", list(puritys_dict.keys()), orientation='vertical', tmp_name = tmp_name)
                # merge_purity_images(f"{output_root_path}/{sample_platform}", f"{output_root_path}/{sample_platform}", list(puritys_dict.keys()), orientation='horizontal', tmp_name = tmp_name)
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

def plot_ratio_t_distribution(somatic_df, output_path, column_name="ratio_t"):
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
    if "HIGH" in column_name:
        bins = np.linspace(0, somatic_df[column_name].max(), 21)  # 0到1之间创建10个区间
    
    # 绘制tp的直方图
    ax1.hist(tp_df[column_name], bins=bins, alpha=0.7, color='blue', edgecolor='black')
    ax1.set_title('True Positives (TP)')
    ax1.set_xlabel(column_name)
    ax1.set_ylabel('Count')
    ax1.grid(True, linestyle='--', alpha=0.7)
    # ax1.set_ylim(0, 2500)  # 固定 y 轴高度
    
    # 绘制fp的直方图
    ax2.hist(fp_df[column_name], bins=bins, alpha=0.7, color='red', edgecolor='black')
    ax2.set_title('False Positives (FP)')
    ax2.set_xlabel(column_name)
    ax2.set_ylabel('Count')
    ax2.grid(True, linestyle='--', alpha=0.7)
    # ax2.set_ylim(0, 2500)  # 固定 y 轴高度
    
    # 添加统计信息
    # tp_mean = tp_df[column_name].mean()
    # tp_median = tp_df[column_name].median()
    # fp_mean = fp_df[column_name].mean()
    # fp_median = fp_df[column_name].median()
    
    # ax1.text(0.05, 0.95, f'Mean: {tp_mean:.3f}\nMedian: {tp_median:.3f}\nCount: {len(tp_df)}', 
    #          transform=ax1.transAxes, verticalalignment='top',
    #          bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    # ax2.text(0.05, 0.95, f'Mean: {fp_mean:.3f}\nMedian: {fp_median:.3f}\nCount: {len(fp_df)}', 
    #          transform=ax2.transAxes, verticalalignment='top',
    #          bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    # 调整布局并保存图片
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()
    
    print(f"已生成ratio_t分布图，保存在 {output_path}")

def plot_af_distribution(test_df, output_path, file_suffix='_ans_vaf'):
    """
    绘制AF值分布统计图表
    
    参数:
        test_df: 包含AF列的DataFrame
        output_path: 输出图片的基础路径
        file_suffix: 文件名后缀，默认为'_ans_vaf'
    """
    if test_df.empty or 'AF' not in test_df.columns:
        logging.warning("DataFrame为空或不包含AF列，无法生成图表")
        return
        
    # 确保输出目录存在
    output_file = f"{output_path}{file_suffix}.png"
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    
    # 处理AF数据
    test_df['AF'] = test_df['AF'].astype(float)
    test_df['AF'] = test_df['AF'].round(2)
    
    # 统计AF值的分布
    af_counts = test_df['AF'].value_counts().sort_index()

    # 创建图表
    plt.figure(figsize=(12, 6))

    # 绘制柱状图
    plt.bar(af_counts.index, af_counts.values, color='skyblue', width=0.008)

    plt.title('AF Distribution for True Positives')
    plt.xlabel('Allele Frequency (AF)')
    plt.ylabel('Count')

    # 设置x轴范围和刻度
    plt.xlim(0, 1)
    x_ticks = np.arange(0, 1.1, 0.05)
    plt.xticks(x_ticks, [f"{tick:.2f}" for tick in x_ticks], rotation=45, ha='right')

    plt.grid(True, axis='y', linestyle='--', alpha=0.7)

    plt.tight_layout()

    # 保存图表
    plt.savefig(output_file, dpi=300)
    plt.close()
    
    logging.info(f"已生成AF值分布统计图表，保存在 {output_file}")

def plot_value_distribution(tp_df, fp_df, output_path, title, value_column):

    """
    绘制 TP 和 FP 的 HIGHV 分布条形图
    
    参数:
        tp_df: 包含 HIGHV 列的 TP DataFrame
        fp_df: 包含 HIGHV 列的 FP DataFrame
        output_path: 输出图片的路径
        title: 图表标题
        value_column: 要统计的列名
    """
    import matplotlib.pyplot as plt
    import numpy as np
    import os
    
    # 确保输出目录存在
    output_dir = os.path.dirname(output_path)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # 获取值的分布
    tp_counts = tp_df[value_column].value_counts().sort_index()
    fp_counts = fp_df[value_column].value_counts().sort_index()
    
    # 获取所有可能的值
    all_values = sorted(set(tp_counts.index) | set(fp_counts.index))
    
    # 创建图表
    plt.figure(figsize=(12, 6))
    
    # 设置条形图的位置
    x = np.arange(len(all_values))
    width = 0.35
    
    # 绘制条形图
    tp_bars = plt.bar(x - width/2, [tp_counts.get(val, 0) for val in all_values], 
                     width, label='True Positives', color='blue', alpha=0.7)
    fp_bars = plt.bar(x + width/2, [fp_counts.get(val, 0) for val in all_values], 
                     width, label='False Positives', color='red', alpha=0.7)
    
    # 添加标题和标签
    plt.title(title)
    plt.xlabel(f'{value_column} Value')
    plt.ylabel('Count')
    
    # 设置 x 轴刻度，每隔 5 个值显示一个标签
    plt.xticks(x[::5], all_values[::5], rotation=45)
    
    plt.legend()
    
    # 添加网格线
    plt.grid(True, axis='y', linestyle='--', alpha=0.7)
    
    # 添加统计信息
    tp_total = len(tp_df)
    fp_total = len(fp_df)
    plt.text(0.02, 0.95, f'TP Total: {tp_total}\nFP Total: {fp_total}', 
             transform=plt.gca().transAxes, verticalalignment='top',
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    # 调整布局并保存图片
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()
    
    print(f"已生成 {value_column} 分布图，保存在 {output_path}")


def plot_ssrs_votes_comparison(vote_df, high_cols, low_cols, precision_calculator,get_result_name, variant_type, calling_vcf_path, ans_vcf_path, ans_bed_path, run_chr, output_path):
    """Generate and save the SSRS votes comparison bar chart."""
    # Calculate total votes
    vote_df["total_high"] = vote_df[high_cols].sum(axis=1)
    vote_df["total_low"] = vote_df[low_cols].sum(axis=1)

    # Filter rows with exactly one high vote
    # vote_df.query("total_high == 1", inplace=True)

    # Get result and filter TP and FP
    vote_ans_df = get_result(
        precision_calculator, 
        get_result_name,
        in_vcf_path=calling_vcf_path, 
        in_vcf_df=vote_df, 
        ans_vcf_path=ans_vcf_path, 
        ans_bed_path=ans_bed_path,
        test_chr=run_chr, 
        check_type=variant_type
    ).query("ans == 'tp' or ans == 'fp'")

    # Separate TP and FP dataframes
    vote_ans_tp_df = vote_ans_df.query("ans == 'tp'")
    vote_ans_fp_df = vote_ans_df.query("ans == 'fp'")

    # Calculate total high votes
    total_high_votes_tp = vote_ans_tp_df["total_high"].sum()
    total_high_votes_fp = vote_ans_fp_df["total_high"].sum()

    # Calculate votes and percentages
    tp_votes = [sum(vote_ans_tp_df[col]) for col in high_cols]
    tp_percentages = [f"{votes/total_high_votes_tp:.2%}" for votes in tp_votes]
    fp_votes = [sum(vote_ans_fp_df[col]) for col in high_cols]
    fp_percentages = [f"{votes/total_high_votes_fp:.2%}" for votes in fp_votes]

    print("TP votes:", *tp_votes)
    print("TP percentages:", *tp_percentages)
    print("FP votes:", *fp_votes)
    print("FP percentages:", *fp_percentages)

    # Calculate fractions for plotting
    tp_frac = [votes/total_high_votes_tp for votes in tp_votes]
    fp_frac = [votes/total_high_votes_fp for votes in fp_votes]

    # Plotting the bar chart
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.ticker as mtick

    fig, ax = plt.subplots(figsize=(10, 6))
    x = np.arange(len(high_cols))
    width = 0.35

    ax.bar(x - width/2, tp_frac, width, label='TP', color='green')
    ax.bar(x + width/2, fp_frac, width, label='FP', color='blue')
    ax.set_xticks(x)
    ax.set_xticklabels(high_cols, rotation=45)
    ax.set_ylabel('Percentage')
    ax.set_title('TP and FP Votes Comparison')
    ax.legend()
    ax.yaxis.set_major_formatter(mtick.PercentFormatter(1.0))

    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()

def test_precision():
    precision_calculator = PrecisionCalculator()
    # calling_vcf_path = "/big8_disk/data/HCC1395/ONT/ClairS_TO_v0_3_0/indel.vcf"
    calling_vcf_path = "/big8_disk/fenne113/clairs-to_output_ssrs/tmp/vcf_output/indel_pileup_nonsomatic_tagging.vcf"
    calling_vcf_path = "/big8_disk/fenne113/clairs-to_output_ssrs/indel.vcf"
    calling_vcf_path = "/bip8_disk/fenne113/indel/tumor/snp_with_tp_tagging_sorted.vcf"

    ans_vcf_path = "/big8_disk/data/HCC1395/SEQC2/high-confidence_sINDEL_in_HC_regions_v1.2.1.vcf"
    ans_bed_path = "/big8_disk/data/HCC1395/SEQC2/High-Confidence_Regions_v1.2.bed"
    variant_type = "indel"
    get_result(precision_calculator, f"test",
                in_vcf_path=calling_vcf_path, in_vcf_df=None, ans_vcf_path=ans_vcf_path, ans_bed_path=ans_bed_path, check_type=variant_type)
    precision_calculator.print_confusion_matrix_metrics()

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
    # test_precision()
    