import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor, as_completed

# 腫瘤純度 ρ 取 0.1 到 1.0 (間隔 0.05)
rho_values = np.arange(0.1, 1.01, 0.05)
# 腫瘤倍性 ψₜ 取 1.6 到 4.8 (間隔 0.1)
psi_values = np.arange(1.6, 4.81, 0.1)

def compute_nA_nB(vaf, rho, psi_t, r=0, gamma=1):
    psi = 2 * (1 - rho) + rho * psi_t
    x = psi * (2 ** (r / gamma))
    nA = (x * (1 - vaf) - (1 - rho)) / rho
    nB = (x * vaf - (1 - rho)) / rho
    return nA, nB

def compute_error(nA, nB, vaf):
    err = (nA - np.round(nA))**2 + (nB - np.round(nB))**2
    
    # weight = 0.05 if np.abs(vaf - 0.5) < 1e-3 else 1.0
    weight = np.where(np.abs(vaf - 0.5) < 1e-3, 0.05, 1.0)  # 修正錯誤

    return weight * err

def compute_percentzero(nA, nB, lengths):
    mask_nA = np.round(nA) == 0  # 判斷 nA 是否接近 0
    mask_nB = np.round(nB) == 0  # 判斷 nB 是否接近 0
    
    numerator = np.sum(mask_nA * lengths) + np.sum(mask_nB * lengths)
    denominator = np.sum(lengths)
    
    return numerator / denominator if denominator != 0 else 0  # 避免除以零

def create_distance_matrix(baf):

    # # 腫瘤純度 ρ 取 0.1 到 1.0 (間隔 0.05)
    # rho_values = np.arange(0.1, 1.01, 0.05)
    # # 腫瘤倍性 ψₜ 取 1.6 到 4.8 (間隔 0.1)
    # psi_values = np.arange(1.6, 4.81, 0.1)

    # d = np.zeros((len(psi_values), len(rho_values)))
    d = np.zeros((len(psi_values), len(rho_values)), dtype=object)

    MINRHO = 0.2
    MINPERCZEROABB = 0.1
    MINPLOIDYSTRICT = 1.7
    MAXPLOIDYSTRICT = 2.3

    for i, psi in enumerate(psi_values):
        for j, rho in enumerate(rho_values):
            # nA = (rho - 1 - (s[:, 1] - 1) * 2 ** (s[:, 0] / gamma) * ((1 - rho) * 2 + rho * psi)) / rho
            # nB = (rho - 1 + s[:, 1] * 2 ** (s[:, 0] / gamma) * ((1 - rho) * 2 + rho * psi)) / rho
            nA, nB = compute_nA_nB(baf, rho, psi)
            err = compute_error(nA, nB, baf)
            ploidy = sum((nA+nB) * baf) / sum(baf)
            # if (ploidy > MINPLOIDYSTRICT and ploidy < MAXPLOIDYSTRICT and rho >= MINRHO and err > MINPERCZEROABB):
            #     d[i, j] = err
            mask = (MINPLOIDYSTRICT < ploidy) & (ploidy < MAXPLOIDYSTRICT) & (rho >= MINRHO)
            err_mask = np.where(err > MINPERCZEROABB, err, np.nan)  # 只保留符合條件的 err，其他設為 NaN

            if mask:
                d[i, j] = err_mask  # 儲存篩選後的 err 陣列

            # percentzero = compute_percentzero(nA, nB, lengths)
            # if np.nansum(nA) < np.nansum(nB):
            #     nMinor = nA
            # else:
            #     nMinor = nB

            # d[i, j] = np.nansum(
            #     np.abs(nMinor - np.maximum(np.round(nMinor), 0)) ** 2 * s[:, 2] * np.where(s[:, 1] == 0.5, 0.05, 1))

    return d

def plot_heatmap(distance_matrix, name):
    """
    Plot heatmap of distance matrix
    
    Args:
        distance_matrix: The matrix containing error values
        rho_values: Array of rho values (tumor purity)
        psi_values: Array of psi values (tumor ploidy)
    """
    # Convert object array to numeric array for plotting
    numeric_matrix = np.zeros((len(psi_values), len(rho_values)))
    
    for i in range(len(psi_values)):
        for j in range(len(rho_values)):
            if distance_matrix[i, j] is not None:
                # If the cell contains an array, take the mean of non-NaN values
                if isinstance(distance_matrix[i, j], np.ndarray):
                    valid_values = distance_matrix[i, j][~np.isnan(distance_matrix[i, j])]
                    numeric_matrix[i, j] = np.mean(valid_values) if len(valid_values) > 0 else np.nan
                else:
                    numeric_matrix[i, j] = distance_matrix[i, j]
            else:
                numeric_matrix[i, j] = np.nan
    
    # Create figure
    plt.figure(figsize=(12, 10))
    
    # Create heatmap
    im = plt.imshow(numeric_matrix, cmap='viridis', aspect='auto', origin='lower')
    
    # Add colorbar
    plt.colorbar(im, label='Error')
    
    # Set ticks and labels
    plt.xticks(np.arange(0, len(rho_values), 2), [f'{rho:.2f}' for rho in rho_values[::2]])
    plt.yticks(np.arange(0, len(psi_values), 2), [f'{psi:.1f}' for psi in psi_values[::2]])
    
    # Add labels and title
    plt.xlabel('tumor purity (ρ)')
    plt.ylabel('tumor ploidy (ψₜ)')
    plt.title('error heatmap')
    
    # Add grid
    plt.grid(False)
    
    plt.savefig(name)
    # Show plot
    plt.tight_layout()
    plt.show()

# 主程式範例
if __name__ == '__main__':
    # 假設一個 vaf 值（BAF），例如 0.65
    vaf_example = np.array([0.26])
    test = create_distance_matrix(vaf_example)
    
    
    # 繪製熱力圖
    plot_heatmap(test, 'test.png')
    
    print(test)

    # 若有 r 與 gamma 資訊，可於此指定，目前 r=0, gamma=1
    # results_df, best = grid_search(vaf_example, r=0, gamma=1)
    # print(results_df)
    # print(best)