from utils import *
from config_parser import get_config_column
import matplotlib.patches as patches
from typing import List, Tuple, Dict, Optional
import pandas as pd
import matplotlib.pyplot as plt
from dataclasses import dataclass
from enum import Enum

class BandType(Enum):
    """染色體帶型類型"""
    GNEG = "gneg"      # 白色
    GPOS25 = "gpos25"  # 淺灰
    GPOS50 = "gpos50"  # 中灰
    GPOS75 = "gpos75"  # 深灰
    GPOS100 = "gpos100" # 黑色
    ACEN = "acen"      # 紅色(著絲粒)
    GVAR = "gvar"      # 變異區域
    STALK = "stalk"    # stalk區域

@dataclass
class Region:
    """基因組區域數據類"""
    chromosome: str
    start_pos: int
    end_pos: int
    
    @classmethod
    def from_string(cls, region_str: str) -> 'Region':
        """從字符串創建Region對象"""
        try:
            if not region_str:
                return cls("all", -1, -1)
                
            if ':' not in region_str:
                return cls(region_str, -1, -1)
                
            chrom, pos = region_str.split(':')
            
            if '-' not in pos:
                start_pos = int(pos.replace(',', ''))
                return cls(chrom, start_pos, start_pos)
                
            start, end = pos.split('-')
            start_pos = int(start.replace(',', ''))
            end_pos = int(end.replace(',', ''))
            
            return cls(chrom, start_pos, end_pos)
            
        except (ValueError, AttributeError):
            raise ValueError("無效的區域格式。預期格式: 'chrX:start-end', 'chrX:pos' 或 'chrX'")

class ChromosomeVisualizer:
    """染色體可視化類"""
    
    def __init__(self, region_str: str = "", padding: int = 10000, output_dir: str = ""):
        """
        初始化染色體可視化器
        
        Args:
            region_str: 基因組區域字符串
            padding: 區域擴展大小
            output_dir: 輸出目錄
        """
        self.region = Region.from_string(region_str)
        self.padding = padding
        self.output_dir = output_dir
        
        # 顏色映射
        self.band_colors = {
            BandType.GNEG.value: "#FFFFFF",    # 白色
            BandType.GVAR.value: "#E0E0E0",    # 變異區域
            BandType.STALK.value: "#C0C0C0",   # stalk區域
            BandType.GPOS25.value: "#C0C0C0",  # 淺灰
            BandType.GPOS50.value: "#808080",  # 中灰
            BandType.GPOS75.value: "#404040",  # 深灰
            BandType.GPOS100.value: "#000000", # 黑色
            BandType.ACEN.value: "#FF0000",    # 紅色(著絲粒)
        }
        
        # 設置要顯示的圖層  # ref, ans, loh, sge, lge
        self.print_bool = [True, False, False, False, False]
        
        # 加載參考數據
        self._load_reference_data()
        
        # 設置繪圖參數
        self._setup_plotting()
        
    def _load_reference_data(self):
        """加載參考數據"""
        ref_bed_file = 'bioinfo_handler/data/hgTables.bed'
        self.ref_bed_data = pd.read_csv(
            ref_bed_file, 
            comment='#', 
            sep='\t', 
            names=["chrom", "chromStart", "chromEnd", "name", "gieStain"]
        )
        
        # 設置染色體列表
        if self.region.chromosome == "all":
            self.chromosomes = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY']
        else:
            self.chromosomes = [self.region.chromosome]
            
        # 設置位置範圍
        if self.region.start_pos == -1:
            self.region.start_pos = 0
            self.region.end_pos = self.ref_bed_data[self.ref_bed_data['chrom'] == 'chr1'].iloc[-1]['chromEnd']
        else:
            self.region.start_pos -= self.padding
            self.region.end_pos += self.padding
            
        self.width = self.region.end_pos - self.region.start_pos
        
    def _setup_plotting(self):
        """設置繪圖參數"""
        self.print_int = sum(self.print_bool)
        
        # 創建子圖
        chrom_height = 1
        n_chroms = len(self.chromosomes)
        height_ratios = [0.15] * self.print_int * n_chroms  # 圖例空間移到最前面
        # height_ratios = [0.05] + [0.15] * self.print_int * n_chroms  # 圖例空間移到最前面
        
        fig, axs = plt.subplots(
            n_chroms * self.print_int,  # 加1為圖例
            # n_chroms * self.print_int + 1,  # 加1為圖例
            1, 
            figsize=(25, n_chroms * self.print_int * chrom_height + 2),  # 增加高度以容納圖例
            gridspec_kw={'height_ratios': height_ratios}
        )
        
        self.chrom_height = chrom_height
        self.axs = [axs] if isinstance(axs, plt.Axes) else axs
        
        plt.subplots_adjust(left=0.10, right=0.99, top=0.95, bottom=0.05, hspace=1, wspace=0)
        
    def _draw_legend(self):
        """繪製統一的圖例"""
        legend_ax = self.axs[0]  # 第一個子圖用於圖例

        # 創建圖例元素
        ref_elements = []
        ref_labels = []
        loh_ans_elements = []
        loh_ans_labels = []
        loh_elements = []
        loh_labels = []

        # 添加參考帶型圖例
        for band_type, color in self.band_colors.items():
            ref_elements.append(patches.Patch(facecolor=color, edgecolor='black', label=band_type))
            ref_labels.append(band_type)

        # 添加LOH分數圖例
        loh_elements.append(patches.Patch(facecolor=self._get_score_color(0), edgecolor='black', label='LOH Score: 0'))
        loh_elements.append(patches.Patch(facecolor=self._get_score_color(0.045), edgecolor='black', label='LOH Score: 0.045'))
        loh_elements.append(patches.Patch(facecolor=self._get_score_color(0.09), edgecolor='black', label='LOH Score: 0.09'))
        loh_labels.extend(['LOH Score: 0', 'LOH Score: 0.045', 'LOH Score: 0.09'])

        # 添加LOH答案圖例
        loh_ans_elements.append(patches.Patch(facecolor=(139/255, 16/255, 136/255), edgecolor='black', label='LOH'))
        loh_ans_elements.append(patches.Patch(facecolor='white', edgecolor='black', label='Non-LOH'))
        loh_ans_labels.extend(['LOH', 'Non-LOH'])

        # 繪製參考帶型圖例（放上面）
        ref_legend = legend_ax.legend(
            handles=ref_elements,
            labels=ref_labels,
            ncol=len(ref_elements),
            loc='upper center',
            bbox_to_anchor=(0.5, 1.05),
            fontsize=12,
            title='Reference Bands',
            title_fontsize=14
        )
        legend_ax.add_artist(ref_legend)

        # 繪製LOH答案圖例（放中間）
        loh_ans_legend = legend_ax.legend(
            handles=loh_ans_elements,
            labels=loh_ans_labels,
            ncol=len(loh_ans_elements),
            loc='center',
            bbox_to_anchor=(0.5, 0.55),
            fontsize=12,
            title='LOH Answer',
            title_fontsize=14
        )
        legend_ax.add_artist(loh_ans_legend)

        # 繪製LOH分數圖例（放下面）
        loh_legend = legend_ax.legend(
            handles=loh_elements,
            labels=loh_labels,
            ncol=len(loh_elements),
            loc='lower center',
            bbox_to_anchor=(0.5, 0.05),
            fontsize=12,
            title='Longphase-to LOH Calling Scores',
            title_fontsize=14
        )

        # 隱藏 legend 軸
        legend_ax.axis('off')

    def _draw_legend_only(self):
        """單獨繪製圖例並保存"""
        # 創建新的圖形
        fig = plt.figure(figsize=(15, 3))
        ax = fig.add_subplot(111)

        # 創建圖例元素
        ref_elements = []
        ref_labels = []
        loh_ans_elements = []
        loh_ans_labels = []
        loh_elements = []
        loh_labels = []

        # 添加參考帶型圖例
        for band_type, color in self.band_colors.items():
            ref_elements.append(patches.Patch(facecolor=color, edgecolor='black', label=band_type))
            ref_labels.append(band_type)

        # 添加LOH分數圖例
        loh_elements.append(patches.Patch(facecolor='white', edgecolor='black', label='Non-LOH'))
        loh_elements.append(patches.Patch(facecolor=self._get_score_color(0.09), edgecolor='black', label='LOH Score: 0.09'))
        loh_elements.append(patches.Patch(facecolor=self._get_score_color(0.045), edgecolor='black', label='LOH Score: 0.045'))
        loh_elements.append(patches.Patch(facecolor=self._get_score_color(0), edgecolor='black', label='LOH Score: 0'))
        loh_labels.extend(['Non-LOH Score: > 0.09', 'LOH Score: 0.09', 'LOH Score: 0.045', 'LOH Score: 0'])

        # 添加LOH答案圖例
        loh_ans_elements.append(patches.Patch(facecolor='white', edgecolor='black', label='Non-LOH'))
        loh_ans_elements.append(patches.Patch(facecolor=(139/255, 16/255, 136/255), edgecolor='black', label='LOH'))
        loh_ans_labels.extend(['Non-LOH', 'LOH'])

        # 繪製參考帶型圖例（放上面）
        ref_legend = ax.legend(
            handles=ref_elements,
            labels=ref_labels,
            ncol=len(ref_elements),
            loc='upper center',
            bbox_to_anchor=(0.5, 0.85),
            fontsize=12,
            title='Reference Bands',
            title_fontsize=14
        )
        ax.add_artist(ref_legend)

        # 繪製LOH答案圖例（放中間）
        loh_ans_legend = ax.legend(
            handles=loh_ans_elements,
            labels=loh_ans_labels,
            ncol=len(loh_ans_elements),
            loc='center',
            bbox_to_anchor=(0.5, 0.5),
            fontsize=12,
            title='LOH Answer',
            title_fontsize=14
        )
        ax.add_artist(loh_ans_legend)

        # 繪製LOH分數圖例（放下面）
        loh_legend = ax.legend(
            handles=loh_elements,
            labels=loh_labels,
            ncol=len(loh_elements),
            loc='lower center',
            bbox_to_anchor=(0.5, 0.15),
            fontsize=12,
            title='Longphase-to LOH Calling Scores',
            title_fontsize=14
        )

        # 隱藏軸
        ax.axis('off')

        # 調整佈局
        plt.tight_layout()

        # 保存圖例
        legend_path = f'{self.output_dir}/legend.png'
        plt.savefig(legend_path, bbox_inches='tight', dpi=300)
        plt.close(fig)
        print(f'已保存圖例到: {legend_path}')

    def draw(self):
        """繪製染色體圖"""
        # 繪製圖例
        # self._draw_legend()
        
        for i, chrom in enumerate(self.chromosomes):
            row_start = i * self.print_int  # 加1因為圖例佔用了第一個位置
            # row_start = i * self.print_int + 1  # 加1因為圖例佔用了第一個位置
            print_row = 0
            
            # 添加染色體標籤
            self._add_chromosome_label(row_start, chrom)
            
            # 繪製參考數據
            if self.print_bool[0]:
                self._plot_reference(row_start, print_row, chrom)
                print_row += 1
            if self.print_bool[1]:
                self._plot_loh_ans(row_start, print_row, chrom)
                print_row += 1
            if self.print_bool[2]:
                self._plot_loh(row_start, print_row, chrom)
                print_row += 1
        
        # 保存主圖
        output_path = f'{self.output_dir}/bed_view.png'
        plt.savefig(output_path, bbox_inches='tight')
        plt.close()
        print(f'已保存圖片到: {output_path}')
        
        # 單獨保存圖例
        self._draw_legend_only()

    def load_loh_ans(self, input_dir: str, ref_bool: bool = True):
        """Read all LOH analysis result images and store them in a dictionary"""
        self.print_bool[1] = True  # Set to True for LOH analysis results
        self._setup_plotting()
        self.ans_images = {}
        # Create regex pattern based on ref_bool
        pattern = r'chr[0-9A-Za-z]+\.png$' if ref_bool else r'chr[0-9A-Za-z]+_loh\.png$'
        for file in os.listdir(input_dir):
            if re.match(pattern, file):
                img_path = os.path.join(input_dir, file)
                key = file.replace('.png', '').replace('_loh', '')  # Remove .png and _loh extensions
                self.ans_images[key] = img_path
        return self
    
    def load_loh(self, input_file: str):
        """Load LOH data from a file"""
        self.print_bool[2] = True  # Set to True for LOH data
        self._setup_plotting()
        self.loh_data = pd.read_csv(input_file, sep='\t', names=['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand', 'thickStart', 'thickEnd', 'itemRgb'])
        return self

    def _add_chromosome_label(self, row_start: int, chrom: str):
        """添加染色體標籤"""
        self.axs[row_start].text(
            -0.02, 0.5, chrom,
            verticalalignment='center',
            horizontalalignment='right',
            fontsize=40,
            fontweight='bold',
            transform=self.axs[row_start].transAxes
        )
        
    def _plot_reference(self, row_start: int, print_row: int, chrom: str):
        """繪製參考數據"""
        chrom_data = self.ref_bed_data[self.ref_bed_data['chrom'] == chrom]
        ax = self.axs[row_start + print_row]
        
        for _, row in chrom_data.iterrows():
            # 檢查是否與目標區域重疊
            band_start = row['chromStart']
            band_end = row['chromEnd']
            
            # 如果完全在目標區域外，跳過
            if band_end < self.region.start_pos or band_start > self.region.end_pos:
                continue
                
            # 計算重疊部分
            draw_start = max(band_start, self.region.start_pos)
            draw_end = min(band_end, self.region.end_pos)
            width = draw_end - draw_start
            
            color = self.band_colors.get(row['gieStain'], "#FFFFFF")
            rect = patches.Rectangle(
                (draw_start, 0),
                width,
                self.chrom_height,
                edgecolor='black',
                facecolor=color
            )
            ax.add_patch(rect)
            
        # 設置軸屬性
        ax.set_xlim(self.region.start_pos, self.region.end_pos)
        ax.get_yaxis().set_visible(False)
        ax.set_yticks([])
        ax.ticklabel_format(useOffset=False, style='plain')
        ax.set_title(f'{chrom} Reference bands', loc='left')
        
    def _plot_loh_ans(self, row_start: int, print_row: int, chrom: str):
        """Plot LOH analysis results with scaled x-axis"""
        ax = self.axs[row_start + print_row]
        img_path = self.ans_images[chrom]
        # Read and plot the image
        img = plt.imread(img_path)
            
        # Calculate scaling factors
        img_width = img.shape[1]
        scale_factor = self.width / img_width
        chrom_min_length = self.ref_bed_data[self.ref_bed_data['chrom'] == chrom]['chromStart'].min()
        chrom_max_length = self.ref_bed_data[self.ref_bed_data['chrom'] == chrom]['chromEnd'].max()

        # 計算圖像中對應的像素範圍
        total_length = chrom_max_length - chrom_min_length
        pixel_per_unit = img_width / total_length
        # 計算需要保留的區域在圖像上的像素範圍
        start_pixel = int((self.region.start_pos - chrom_min_length) * pixel_per_unit)
        end_pixel = int((self.region.end_pos - chrom_min_length) * pixel_per_unit)
        # 切割圖像，只保留指定的區域
        cropped_img = img[:, start_pixel:end_pixel, :]

        start_length = max(self.region.start_pos, chrom_min_length)
        end_length = min(self.region.end_pos, chrom_max_length)
        # Plot the image with proper scaling
        ax.imshow(cropped_img,
                extent=[start_length, end_length, 0, self.chrom_height],
                aspect='auto')# 繪製切割後的圖像，並保持縮放
        # Set axis properties
        ax.set_xlim(self.region.start_pos, self.region.end_pos)
        ax.get_yaxis().set_visible(False)
        ax.set_yticks([])
        ax.ticklabel_format(useOffset=False, style='plain')
        ax.set_title(f'{chrom} LOH Answer', loc='left')

    def _get_score_color(self, score: float) -> str:
        """根據score值返回藍綠漸變色
        
        Args:
            score: 分數值
            
        Returns:
            str: 十六進制顏色代碼
        """
        # 將score歸一化到0-1範圍
        normalized_score = min(max(score, 0), 0.09) / 0.09  # 將0-0.09歸一化到0-1範圍
        r_start, g_start, b_start = 128, 128, 255    # 淺藍
        r_end, g_end, b_end = 128, 255, 128          # 淺綠
        # 使用淺藍綠漸變色
        r = int(r_start + (r_end - r_start) * normalized_score)
        g = int(g_start + (g_end - g_start) * normalized_score)
        b = int(b_start + (b_end - b_start) * normalized_score)
        return f'#{r:02x}{g:02x}{b:02x}'

    def _plot_loh(self, row_start: int, print_row: int, chrom: str):
        """繪製LOH數據，使用藍綠漸變色表示score值"""
        df = self.loh_data[self.loh_data['chrom'] == chrom]
        ax = self.axs[row_start + print_row]
        
        for _, row in df.iterrows():
            # 檢查是否與目標區域重疊
            band_start = row['chromStart']
            band_end = row['chromEnd']
            
            # 如果完全在目標區域外，跳過
            if band_end < self.region.start_pos or band_start > self.region.end_pos:
                continue
                
            # 計算重疊部分
            draw_start = max(band_start, self.region.start_pos)
            draw_end = min(band_end, self.region.end_pos)
            width = draw_end - draw_start
            
            # 根據score值獲取顏色
            color = self._get_score_color(float(row['score']))
            
            # 繪製矩形
            rect = patches.Rectangle(
                (draw_start, 0),
                width,
                self.chrom_height,
                edgecolor='black',
                facecolor=color
            )
            ax.add_patch(rect)
            
        # 設置軸屬性
        ax.set_xlim(self.region.start_pos, self.region.end_pos)
        ax.get_yaxis().set_visible(False)
        ax.set_yticks([])
        ax.ticklabel_format(useOffset=False, style='plain')
        ax.set_title(f'{chrom} Longphase-to LOH Calling', loc='left')

def main():
    """主函數"""
    region_str = 'chr1:5000000-50000000'
    # region_str = 'chr1'
    region_str = ''
    padding = 1000
    root_path = get_config_column('run_root_path')
    plot_folder_name = get_config_column('output_plot_folder')
    folder_name = get_config_column('output_folder')
    output_dir = f"{root_path}/{plot_folder_name}"
    loh_ans_dir = f"{root_path}/{plot_folder_name}/transformed_images2"
    loh_file = f"{root_path}/{folder_name}/HCC1395_ONT/t50_n00/output_LOH.bed"
    
    visualizer = ChromosomeVisualizer(region_str, padding, output_dir)
    visualizer.load_loh_ans(loh_ans_dir, ref_bool=False).load_loh(loh_file)
    # visualizer.load_loh_ans(loh_ans_dir, ref_bool=True).load_loh(loh_file)
    visualizer.draw()

if __name__ == "__main__":
    main()



