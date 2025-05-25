import os
import re
import sys
import glob
import socket
import random
import logging
import subprocess
import numpy as np
import pandas as pd
from PIL import Image
from tqdm import tqdm
import matplotlib.pyplot as plt
from multiprocessing import Process, Pool
from matplotlib.lines import Line2D # Import Line2D for legend handles
import matplotlib.ticker as mtick
from io import StringIO

rainbow_color_map = {
    # "Red": "#e41a1c",       # 紅色
    # "Orange": "#ff7f00",    # 橙色
    # "Yellow": "#ffff33",    # 黃色
    # "Green": "#4daf4a",     # 綠色
    # "Blue": "#377eb8",      # 藍色
    # "Indigo": "#6a3d9a",    # 靛色（偏深紫）
    # "Violet": "#984ea3",    # 紫色
    # "Pink": "#f781bf",      # 粉紅
    # "Brown": "#a65628",     # 棕色
    # "Cyan": "#00ced1",      # 青色

    "Red":    "#d62728",  # 紅色 (Red - 沿用您範例中的紅色)
    "Orange": "#ff7f0e",  # 橙色 (Orange - 沿用您範例中的橘色)
    "Yellow": "#FFD700",  # 黃色 (Yellow - 使用金色 Gold，比純黃 #FFFF00 在白底上更易讀)
    "Green":  "#2ca02c",  # 綠色 (Green - 沿用您範例中的綠色)
    "Blue":   "#1f77b4",  # 藍色 (Blue - 沿用您範例中的藍色)
    "Indigo": "#4B0082",  # 靛色 (Indigo - 傳統的靛藍色)
    "Violet": "#9467bd",  # 紫色 (Violet/Purple - 沿用您範例中的紫色)
    "Pink":   "#e377c2",  # 粉紅色 (Pink - 沿用您範例中的粉紅色)
    "Brown":  "#a65628",  # 棕色 (Brown - 沿用您範例中的棕色)
    "Cyan":   "#00ced1",  # 青色 (Cyan - 沿用您範例中的青色)
    "Gray":   "#808080",  # 灰色 (Gray - 沿用您範例中的灰色)
    "White":  "#ffffff",  # 白色 (White - 沿用您範例中的白色)
    "Black":  "#000000",  # 黑色 (Black - 沿用您範例中的黑色)
}


sample_color_map = {

    "H2009_ONT": "#1f77b4",    # 藍色
    "HCC1954_ONT": "#ff7f0e",  # 橘色
    "HCC1937_ONT": "#2ca02c",  # 綠色
    "COLO829_ONT_PAO": "#d62728",  # 紅色
    "H1437_ONT": "#9467bd",    # 紫色
    "HCC1395_ONT": "#8c564b",   # 棕色
    "HCC1395_ONT_Dorado": "#e377c2"   # 粉紅色
}

softwares = [
    'Longphase_TO_v0_0_1', # 0
    'DeepSomatic_TO_v1_8_0', # 1
    'ClairS_TO_v0_3_0', # 2
    'Longphase_TO_v0_0_1_deepVCF', # 3
    'A', # 4
    'ClairS_TO_ss_v0_3_0', # 5
    'Longphase_TO_ss_v0_0_1', # 7
    'B', # 6
    'ClairS_TO_v0_3_0_pileup_nonsomatic', # 8
    'Longphase_TO_ssrs_v0_0_1', # 9
    'ClairS_TO_ssrs_v0_3_0', # 10
    # 'Longphase_TO_ss_v0_0_1', # 9
    # 'ClairS_TO_ss_v0_3_0', # 10
]

def get_random_color():
    """
    Generate a random color in hex format.
    Returns a random color from the rainbow_color_map.
    """
    return "#{:06x}".format(random.randint(0, 0xFFFFFF))

# Define colors for different software
software_colors = {
    softwares[9]: rainbow_color_map["Orange"],
    softwares[0]: rainbow_color_map["Orange"],
    softwares[1]: rainbow_color_map["Violet"],
    softwares[2]: rainbow_color_map["Blue"],
    softwares[10]: rainbow_color_map["Blue"],
    softwares[3]: rainbow_color_map["Green"],
    softwares[4]: rainbow_color_map["Indigo"],
    softwares[5]: rainbow_color_map["Red"],
    softwares[6]: rainbow_color_map["Cyan"],
    softwares[7]: rainbow_color_map["Brown"],
    softwares[8]: rainbow_color_map["Pink"],
    # softwares[8]: rainbow_color_map["White"],
    'default': get_random_color()
}
markers = ['o', 's', '^', 'D', 'v', '<', '>', 'p', '*', 'h', 'H', '+', 'x', 'd', '|', '_']
software_markers = {
    softwares[0]: markers[0],
    softwares[1]: markers[1],
    softwares[2]: markers[2],
    'default': markers[5]
}

line_styles = {
    softwares[0]: '-',
    softwares[1]: '--',
    softwares[2]: '-.',
    'default': ':'
}
software_rank = {
    softwares[0]: 0,
    softwares[1]: 2,
    softwares[2]: 1,
    'default': 5
}

class AnsiColors:
    """Encapsulate ANSI escape codes for colored output."""
    RESET = "\033[0m"  # Reset to default color
    RED = "\033[31m"
    GREEN = "\033[32m"
    YELLOW = "\033[33m"
    BLUE = "\033[34m"
    PURPLE = "\033[35m"
    CYAN = "\033[36m"
    WHITE = "\033[37m"


def print_all(context):
    with pd.option_context("display.max_rows", None, "display.max_columns", None):
        print(context)


def chr_to_num(chr_name):
    """
    将染色体名称转换为数字。
    例如，将 'chr1' 转换为 1。
    """
    try:
        return int(chr_name.replace("chr", ""))
    except (ValueError, AttributeError):
        return -1  # 或者其他适当的默认值
def split_version_string(version_string):
    """
    Define a function to split the version string
    """
    # Use regex to match the pattern v<number>_<number>_<number>
    match = re.search(r'v(\d+_\d+_\d+)', version_string)
    if match:
        # Split the string at the last underscore
        base_name = version_string[:match.end()]
        suffix = version_string[match.end():]
        return base_name, suffix
    return version_string, None  # Return original if no match


def calculate_purity(name):
    match = re.findall(r't(\d+)_n(\d+)', name)
    if not match:
        return 0
    t, n = match[0]
    purity = round(int(t) / (int(t) + int(n)), 1)
    return purity


def decompress_gzip(file_path):
    subprocess.run(f"gzip -kd {file_path}", shell=True, check=True)


def merge_purity_images(output_path, input_path, input_list, orientation='horizontal', tmp_name = "_ans_vaf"):
    # Collect all purity image paths
    image_paths = [f"{input_path}/{input}{tmp_name}.png" for input in input_list]
    print(image_paths)
    # Open images and store them in a list
    images = [Image.open(image_path)
              for image_path in image_paths if os.path.exists(image_path)]

    if orientation == 'horizontal':
        total_width = sum(image.width for image in images)
        max_height = max(image.height for image in images)

        merged_image = Image.new('RGB', (total_width, max_height))

        current_x = 0
        for image in images:
            merged_image.paste(image, (current_x, 0))
            current_x += image.width

    elif orientation == 'vertical':
        max_width = max(image.width for image in images)
        total_height = sum(image.height for image in images)

        merged_image = Image.new('RGB', (max_width, total_height))

        current_y = 0
        for image in images:
            paste_x = max_width - image.width
            merged_image.paste(image, (paste_x, current_y))
            current_y += image.height

    # Save the merged image
    merged_image_path = f"{output_path}{tmp_name}.png"
    merged_image.save(merged_image_path)
    logging.info(f"Merged purity images saved at: {merged_image_path}")


def estimate_purity(in_df, output_path, name):
    """
    Function to estimate the purity of a sample.
    """
    # in_df = in_df[(in_df['CHROM'] != 'chrX') & (in_df['CHROM'] != 'chrY')]
    # in_df['ratio'] = pd.to_numeric(in_df['ratio'], errors='coerce')
    # in_df = in_df.dropna(subset=['ratio'])  # 移除 NaN
    # print(in_df.query("ans == 'tp' and ratio > 0.98"))
    # print(in_df.query("ans == 'fp' and ratio > 0.98"))
    
    
    # test_df = in_df.query("ans == 'tp'").copy()
    # in_df.drop(columns=['ans'], inplace=True)
    test_df = in_df.copy()
    boxplot_data = test_df['ratio'].describe()

    # Save the boxplot data to the output path
    boxplot_data.to_csv(output_path.replace(".png", "_1.csv"), index=True)

    # bins = np.linspace(0, 1, 101)
    # ratio_bins = pd.cut(in_df['ratio'], bins=bins, include_lowest=True)
    # ratio_counts = ratio_bins.value_counts().sort_index()
    # ratio_counts.to_csv(output_path.replace(".png", "_2.csv"), index=True)

    fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(8, 10))  # 2行1列的子圖

    if 'ans' in in_df.columns:
        # Separate the data into true positives (tp) and false positives (fp)
        tp_data = in_df[in_df['ans'] == 'tp']['ratio']
        fp_data = in_df[in_df['ans'] == 'fp']['ratio']
        # Create a histogram with two colors
        ax[0].hist(tp_data, bins=50, edgecolor='black', alpha=0.7, label='True Positives', color='blue')
        ax[0].hist(fp_data, bins=50, edgecolor='black', alpha=0.7, label='False Positives', color='red')
        ax[0].legend()  # Add legend to differentiate tp and fp
    else:
        # Histogram for ratio (上方)
        ax[0].hist(in_df['ratio'], bins=50, edgecolor='black', alpha=0.7)
    ax[0].set_xlabel('Ratio')
    ax[0].set_ylabel('Count')
    ax[0].set_title(f'Distribution of Ratio {name}')
    ax[0].grid(axis='y', linestyle='--', alpha=0.7)
    # Boxplot for ratio (下方，直立式)
    # ax[1].boxplot(test_df['ratio'], vert=True)  # 改成 vert=True
    # ax[1].set_ylabel('Ratio')  # 調整標籤
    # ax[1].set_title('Boxplot of Ratio')
    # Violin plot for ratio (下方)
    ax[1].violinplot(test_df['ratio'], vert=True)
    ax[1].set_ylabel('Ratio')
    ax[1].set_title('Violin Plot of Ratio')
    ax[1].grid(axis='y', linestyle='--', alpha=0.7)

    plt.tight_layout()  # 調整間距
    plt.savefig(f"{output_path}")  # 儲存圖片
    plt.close()

