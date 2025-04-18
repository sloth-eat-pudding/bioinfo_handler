import numpy as np
import pandas as pd
from itertools import product
from concurrent.futures import ProcessPoolExecutor
import multiprocessing

# 參數範圍
rho_values = np.arange(0.1, 1.01, 0.05)  # 腫瘤純度 0.1 ~ 1.0，間隔 0.05
psi_values = np.arange(1.6, 4.81, 0.1)   # 腫瘤倍性 1.6n ~ 4.8n，間隔 0.1

def calculate_copy_numbers(vaf, rho, psi):
    """
    根據 VAF 計算等位基因拷貝數。
    """
    n_A = (psi / rho) * ((2 - 2*rho) * vaf + rho) / 2
    n_B = (psi / rho) * ((2 - 2*rho) * (1 - vaf) + rho) / 2
    return n_A, n_B

def compute_error(args):
    """
    計算所有 SNP 的誤差，誤差 = (計算值 - 最接近整數值)^2
    """
    rho, psi, df = args
    total_error = 0
    for vaf in df['VAF']:
        n_A, n_B = calculate_copy_numbers(vaf, rho, psi)
        error = (n_A - round(n_A))**2 + (n_B - round(n_B))**2
        total_error += error
    return rho, psi, total_error

def grid_search_best_params(df):
    """
    執行網格搜索以尋找誤差最小的 (rho, psi)，並行處理計算加速
    """
    best_rho, best_psi, min_error = None, None, float('inf')
    param_list = [(rho, psi, df) for rho, psi in product(rho_values, psi_values) if rho * psi != 0]
    
    with ProcessPoolExecutor(max_workers=multiprocessing.cpu_count()) as executor:
        results = executor.map(compute_error, param_list)
        
        for rho, psi, error in results:
            if error < min_error:
                min_error = error
                best_rho, best_psi = rho, psi
    
    return best_rho, best_psi, min_error