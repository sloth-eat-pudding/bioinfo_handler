import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor, as_completed

rho_values = np.arange(0.1, 1.01, 0.05)
psi_values = np.arange(1.6, 4.81, 0.1)

def compute_nA_nB(vaf, rho, psi_t, r=0, gamma=1):
    psi = 2 * (1 - rho) + rho * psi_t
    x = psi * (2 ** (r / gamma))
    nA = (x * (1 - vaf) - (1 - rho)) / rho
    nB = (x * vaf - (1 - rho)) / rho
    return nA, nB

def compute_error(vaf, rho, psi_t, r=0, gamma=1):
    nA, nB = compute_nA_nB(vaf, rho, psi_t, r, gamma)
    err = (nA - np.round(nA))**2 + (nB - np.round(nB))**2
    
    weight = 0.05 if np.abs(vaf - 0.5) < 1e-3 else 1.0
    return weight * err, nA, nB

def compute_grid_entry(params):
    vaf, rho, psi_t, r, gamma = params
    err, nA, nB = compute_error(vaf, rho, psi_t, r, gamma)
    return {'rho': rho, 'psi_t': psi_t, 'error': err, 'nA': nA, 'nB': nB}

def grid_search_parallel(vaf, r=0, gamma=1, max_workers=None):
    params = [(vaf, rho, psi_t, r, gamma) for rho in rho_values for psi_t in psi_values]

    if max_workers is None:
        max_workers = max(1, mp.cpu_count() - 1)  # 保留 1 核心，避免完全卡死

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        results = list(executor.map(compute_grid_entry, params, chunksize=len(params) // max_workers))

    results_df = pd.DataFrame(results)
    best_row = results_df.loc[results_df['error'].idxmin()]
    return results_df, best_row

def process_sample(vaf, r=0, gamma=1, max_workers=None):
    _, best_sol = grid_search_parallel(vaf, r, gamma, max_workers)
    return best_sol

def gomain(sample_vaf_df, purity_image_path):
    best_params = []
    max_workers = 30
    # 改用 ProcessPoolExecutor 避免 daemon 問題
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(process_sample, vaf, max_workers=max_workers): vaf for vaf in sample_vaf_df['VAF']}
        for future in as_completed(futures):
            best_params.append(future.result())

    best_params_df = pd.DataFrame(best_params)
    plot_purity_distribution(best_params_df, purity_image_path)

def plot_purity_distribution(df, output_path):
    plt.figure(figsize=(8, 6))
    plt.hist(df['rho'], bins=20, edgecolor='k', alpha=0.7)
    plt.xlabel("Tumor Purity (ρ)")
    plt.ylabel("Frequency")
    plt.title("Tumor Purity Distribution")
    plt.savefig(output_path)
    plt.close()


# 主程式範例
if __name__ == '__main__':
    # 設定網格參數：
    # 腫瘤純度 ρ 取 0.1 到 1.0 (間隔 0.05)
    rho_values = np.arange(0.1, 1.01, 0.05)
    # 腫瘤倍性 ψₜ 取 1.6 到 4.8 (間隔 0.1)
    psi_values = np.arange(1.6, 4.81, 0.1)
    
    # 假設一個 vaf 值（BAF），例如 0.65
    vaf_example = 0.65
    # 若有 r 與 gamma 資訊，可於此指定，目前 r=0, gamma=1
    results_df, best = grid_search_parallel(vaf_example, r=0, gamma=1)
    print(results_df)
    print(best)
    
    # print("對 vaf = {}，最佳參數為：".format(vaf_example))
    # print("  ρ = {:.3f}, ψₜ = {:.3f}".format(best['rho'], best['psi_t']))
    # print("  推導 n_A = {:.3f}, n_B = {:.3f}".format(best['nA'], best['nB']))
    # print("  誤差值 = {:.3f}".format(best['error']))
    
    # # 若有一個 dataframe，其中包含多筆 vaf 值，則可逐筆計算最佳解，並最後看 ρ 分布
    # # 例如，建立一個範例 dataframe
    # sample_vaf_df = pd.DataFrame({
    #     'vaf': [0.65, 0.52, 0.48, 0.70, 0.30, 0.5]  # 注意：vaf 約等於 0.5 的會給予較小權重
    # })
    
    # best_params = []
    # for idx, row in sample_vaf_df.iterrows():
    #     _, best_sol = grid_search(row['vaf'], rho_values, psi_values, r=0, gamma=1)
    #     best_params.append(best_sol)
    
    # best_params_df = pd.DataFrame(best_params)
    # print("\n多筆 vaf 的最佳參數：")
    # print(best_params_df)
    
    # # 繪製最佳 ρ 值的分布圖
    # plt.figure(figsize=(6,4))
    # plt.hist(best_params_df['rho'], bins=len(rho_values), edgecolor='black')
    # plt.xlabel('ρ (腫瘤純度)')
    # plt.ylabel('頻數')
    # plt.title('最佳 ρ 分布')
    # plt.show()
