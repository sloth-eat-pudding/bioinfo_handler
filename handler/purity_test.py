#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Purity simulation module for bioinformatics analysis.
This module simulates read depth and variant allele frequency (VAF) based on tumor purity.
"""

import numpy as np
import matplotlib.pyplot as plt
from typing import Tuple, Dict, List, Optional


class PuritySimulator:
    """
    A class to simulate read depth and VAF based on tumor purity.
    """
    
    def __init__(self, depth: int = 50, purity: float = 0.3):
        """
        Initialize the PuritySimulator.
        
        Args:
            depth: Total sequencing depth
            purity: Tumor purity (0.0 to 1.0)
        """
        self.depth = depth
        self.purity = purity
        self.tumor_reads = int(depth * purity)
        self.normal_reads = int(depth * (1 - purity))
        
    def simulate_normal_haplotype(self, na_ratio: float = 0.5, nb_ratio: float = 0.5) -> Dict[str, int]:
        """
        Simulate normal haplotype reads.
        
        Args:
            na_ratio: Ratio of haplotype A in normal sample
            nb_ratio: Ratio of haplotype B in normal sample
        
        Returns:
            Dictionary with haplotype read counts
        """
        if abs(na_ratio + nb_ratio - 1.0) > 0.001:
            raise ValueError("Haplotype ratios must sum to 1.0")
            
        hp1_reads = int(self.normal_reads * na_ratio)
        hp2_reads = int(self.normal_reads * nb_ratio)
        
        # Adjust for rounding errors
        hp1_reads = self.normal_reads - hp2_reads
        
        return {
            "hp1": hp1_reads,
            "hp2": hp2_reads
        }
    
    def simulate_tumor_haplotype(self, na_ratio: float, nb_ratio: float) -> Dict[str, int]:
        """
        Simulate tumor haplotype reads.
        
        Args:
            na_ratio: Ratio of haplotype A in tumor sample
            nb_ratio: Ratio of haplotype B in tumor sample
        
        Returns:
            Dictionary with haplotype read counts
        """
        if abs(na_ratio + nb_ratio - 1.0) > 0.001:
            raise ValueError("Haplotype ratios must sum to 1.0")
            
        hp1_reads = int(self.tumor_reads * na_ratio)
        hp2_reads = int(self.tumor_reads * nb_ratio)
        
        # Adjust for rounding errors
        hp1_reads = self.tumor_reads - hp2_reads
        
        return {
            "hp1": hp1_reads,
            "hp2": hp2_reads
        }


def simulate_all_purities(depth: int = 50, vaf: float = 0.5, alt_hp: str = 'hp1', 
                         TnA: float = 1/3, TnB: float = 2/3, 
                         purities: List[float] = None):
    """
    模拟不同纯度下的读取分配，并创建直条图。
    
    Args:
        depth: 总测序深度
        vaf: 变异等位基因频率
        alt_hp: 变异所在的单倍型 ('hp1' 或 'hp2')
        TnA: 肿瘤样本中单倍型A的比例
        TnB: 肿瘤样本中单倍型B的比例
        purities: 要模拟的纯度列表，默认为[0.1, 0.3, 0.5, 0.7, 0.9]
    """
    if purities is None:
        purities = np.linspace(0.1, 1.0, 10)  # 最穩

    # 计算变异读取数
    alt_reads = int(depth * vaf)
    
    # 创建图形
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # 颜色定义
    colors = {
        'normal_hp1': 'lightgreen',
        'normal_hp2': 'lightblue',
        'tumor_hp1': 'darkgreen',
        'tumor_hp2': 'darkblue',
        'alt_hp1': 'orange',
        'alt_hp2': 'purple'
    }
    
    # 存储每个纯度的总深度
    total_depths = []
    # 需要存各區段資料來堆疊
    normal_hp1_list = []
    normal_hp2_list = []
    tumor_hp1_list = []
    tumor_hp2_list = []
    alt_hp1_list = []
    alt_hp2_list = []
    fake_hp_list = []
    for purity in purities:
        simulator = PuritySimulator(depth=depth, purity=purity)
        
        # 模拟正常样本单倍型
        normal_hp = simulator.simulate_normal_haplotype(0.5, 0.5)
        
        # 模拟肿瘤样本单倍型
        tumor_hp = simulator.simulate_tumor_haplotype(round(TnA, 3), round(TnB, 3))
        
        # 保存各區段
        fake_hp_list.append(tumor_hp[alt_hp] < alt_reads)  # 修正此行
        if fake_hp_list[-1]:  # Check if the last entry in fake_hp_list is True
            tumor_hp['hp1'] = 0  # Set tumor_hp1 to 0
            tumor_hp['hp2'] = 0  # Set tumor_hp2 to 0
            normal_hp['hp1'] = 0  # Set normal_hp1 to 0
            normal_hp['hp2'] = 0  # Set normal_hp2 to 0
            alt_hp1_list.append(0)
            alt_hp2_list.append(0)
        else:
            if alt_hp == 'hp1':
                alt_hp1_list.append(alt_reads)
                alt_hp2_list.append(0)
                tumor_hp['hp1'] = tumor_hp['hp1'] - alt_reads
            else:
                alt_hp1_list.append(0)
                alt_hp2_list.append(alt_reads)
                tumor_hp['hp2'] = tumor_hp['hp2'] - alt_reads
        normal_hp1_list.append(normal_hp['hp1'])
        normal_hp2_list.append(normal_hp['hp2'])
        tumor_hp1_list.append(tumor_hp['hp1'])
        tumor_hp2_list.append(tumor_hp['hp2'])
        print(f"purity: {purity}")
        # print(f"fake_hp: {fake_hp_list[-1]}")  # 只打印当前的fake_hp
        print(f"normal_hp1_list: {normal_hp1_list}")
        print(f"normal_hp2_list: {normal_hp2_list}")
        print(f"tumor_hp1_list: {tumor_hp1_list}")
        print(f"tumor_hp2_list: {tumor_hp2_list}")
        print(f"alt_hp1_list: {alt_hp1_list}")
        print(f"alt_hp2_list: {alt_hp2_list}")

    # 開始畫圖
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # 堆疊第一層 (normal_hp1)
    ax.bar(purities, normal_hp1_list, 
        width=0.08, align='center', 
        label='Normal HP1', alpha=1, color=colors['normal_hp1'])
    ax.bar(purities, normal_hp2_list,
           bottom=normal_hp1_list, 
           width=0.08, align='center', 
           label='Normal HP2', alpha=1, color=colors['normal_hp2'])
    normal_total = np.array(normal_hp1_list) + np.array(normal_hp2_list)
    ax.bar(purities, tumor_hp1_list,
           bottom=normal_total,
           width=0.08, align='center',
           label='Tumor HP1', alpha=0.7, color=colors['tumor_hp1'])
    ax.bar(purities, tumor_hp2_list,
           bottom=normal_total + tumor_hp1_list,
           width=0.08, align='center',
           label='Tumor HP2', alpha=0.7, color=colors['tumor_hp2'])
    tumor_total = np.array(tumor_hp1_list) + np.array(tumor_hp2_list)+np.array(normal_total)
    ax.bar(purities, alt_hp1_list,
           bottom=tumor_total,
           width=0.08, align='center',
           label='Alt HP1', alpha=0.7, color=colors['alt_hp1'])
    ax.bar(purities, alt_hp2_list,
           bottom=tumor_total + alt_hp1_list,
           width=0.08, align='center',
           label='Alt HP2', alpha=0.7, color=colors['alt_hp2'])

    ax.set_title('Stacked Bar Chart of Haplotype Reads vs Tumor Purity')
    ax.set_xlabel('Tumor Purity')
    ax.set_ylabel('Reads Count')
    ax.set_xticks(purities)
    ax.set_xticklabels([f"{p:.1f}" for p in purities])

    ax.legend()
    
    ax.grid(axis='y', linestyle='--', alpha=1)
    plt.tight_layout()
    plt.show()
    plt.savefig(f'purity_simulation.png')
    
    return fig


def main():
    """
    Main function to demonstrate the PuritySimulator.
    """
    # 固定變因
    depth = 50
    TnA = 2
    TnB = 1
    TnA_ratio = TnA/(TnA+TnB)
    TnB_ratio = TnB/(TnA+TnB)
    vaf = 0.3
    # 測試參數
    purity = 0.3
    alt_hp = 'hp1'  # 假设变异在HP1上
    
    print(f"模拟参数:")
    print(f"总深度: {depth}x")
    print(f"肿瘤纯度: {purity}")
    
    simulator = PuritySimulator(depth=depth, purity=purity)
    
    # 模拟正常样本单倍型
    normal_hp = simulator.simulate_normal_haplotype(0.5, 0.5)
    print(f"\n正常样本单倍型 (1:1):")
    print(f"HP1: {normal_hp['hp1']} reads")
    print(f"HP2: {normal_hp['hp2']} reads")
    
    tumor_hp = simulator.simulate_tumor_haplotype(round(TnA_ratio, 3), round(TnB_ratio, 3))
    print(f"\n肿瘤样本单倍型 ({TnA_ratio:.3f}:{TnB_ratio:.3f}):")
    print(f"HP1: {tumor_hp['hp1']} reads")
    print(f"HP2: {tumor_hp['hp2']} reads")
    
    alt_reads = depth*vaf
    fake_hp = tumor_hp[alt_hp] < alt_reads
    print(f"fake_hp: {fake_hp}")
    
    # 模拟所有纯度并创建直条累积图
    simulate_all_purities(depth=depth, vaf=vaf, alt_hp=alt_hp, TnA=TnA_ratio, TnB=TnB_ratio)


if __name__ == "__main__":
    main()
