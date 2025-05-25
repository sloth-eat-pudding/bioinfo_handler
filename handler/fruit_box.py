import time
import random  # 添加random模塊
import concurrent.futures  # 添加並行處理模塊
from concurrent.futures import ProcessPoolExecutor, as_completed
from numpy import size
from tqdm import tqdm

def calculate_box_sum(grid, r1, c1, r2, c2):
    """
    計算矩形區域內數字的總和。
    r1, c1: 左上角座標 (row, column)
    r2, c2: 右下角座標 (row, column)
    """
    total_sum = 0
    # 確保座標在網格範圍內
    if not (0 <= r1 <= r2 < len(grid) and 0 <= c1 <= c2 < len(grid[0])):
        return -1 # 無效的框選範圍

    for r in range(r1, r2 + 1):
        for c in range(c1, c2 + 1):
            total_sum += grid[r][c]
    return total_sum

# def find_valid_boxes(grid):
#     """
#     找出所有數字總和為 10 的矩形框選。
#     返回一個列表，每個元素是一個字典，包含框選的座標、分數和包含的單元格座標。
#     """
#     valid_boxes = []
#     rows = len(grid)
#     cols = len(grid[0])

#     # 遍歷所有可能的左上角座標 (r1, c1)
#     for r1 in range(rows):
#         for c1 in range(cols):
#             # 遍歷所有可能的右下角座標 (r2, c2)，確保 r2 >= r1 且 c2 >= c1
#             for r2 in range(r1, rows):
#                 for c2 in range(c1, cols):
#                     box_sum = calculate_box_sum(grid, r1, c1, r2, c2)

#                     # 如果總和為 10，則為一個有效框選
#                     if box_sum == 10:
#                         # 計算框選的分數 (包含的單元格數量)
#                         score = (r2 - r1 + 1) * (c2 - c1 + 1)
#                         # 記錄框選包含的所有單元格座標，方便後續檢查是否重疊
#                         cells = [(r, c) for r in range(r1, r2 + 1) for c in range(c1, c2 + 1)]
#                         valid_boxes.append({
#                             'coords': (r1, c1, r2, c2),
#                             'score': score,
#                             'cells': cells
#                         })
#     return valid_boxes

def count_nonzero(grid, r1, c1, r2, c2):
    """統計矩形內值≠0的格子數 (作為分數)。"""
    return sum(1 for r in range(r1, r2 + 1)
                 for c in range(c1, c2 + 1)
                 if grid[r][c] != 0)
def find_valid_boxes(grid):
    """列出所有『和=10』的矩形，附帶「目前有效格子數」當作分數。
    先找2個數字的框，找不到才找3個數字的框，以此類推。"""
    rows, cols = len(grid), len(grid[0])
    boxes = []
    
    # 先找2個數字的框
    for r1 in range(rows):
        for c1 in range(cols):
            for r2 in range(r1, rows):
                for c2 in range(c1, cols):
                    # 計算框內格子數
                    cell_count = (r2 - r1 + 1) * (c2 - c1 + 1)
                    # 只考慮2個數字的框
                    if cell_count == 2:
                        if calculate_box_sum(grid, r1, c1, r2, c2) == 10:
                            score = count_nonzero(grid, r1, c1, r2, c2)
                            if score:  # 至少有 1 格還不是 0
                                cells = [(r, c) for r in range(r1, r2 + 1)
                                                 for c in range(c1, c2 + 1)]
                                boxes.append({
                                    'coords': (r1, c1, r2, c2),
                                    'score': score,
                                    'cells': cells
                                })
    
    # 如果找到2個數字的框，直接返回
    if boxes:
        return boxes
    
    # 如果沒找到2個數字的框，找3個數字的框
    for r1 in range(rows):
        for c1 in range(cols):
            for r2 in range(r1, rows):
                for c2 in range(c1, cols):
                    # 計算框內格子數
                    cell_count = (r2 - r1 + 1) * (c2 - c1 + 1)
                    # 只考慮3個數字的框
                    if cell_count == 3:
                        if calculate_box_sum(grid, r1, c1, r2, c2) == 10:
                            score = count_nonzero(grid, r1, c1, r2, c2)
                            if score:  # 至少有 1 格還不是 0
                                cells = [(r, c) for r in range(r1, r2 + 1)
                                                 for c in range(c1, c2 + 1)]
                                boxes.append({
                                    'coords': (r1, c1, r2, c2),
                                    'score': score,
                                    'cells': cells
                                })
    
    # 如果找到3個數字的框，直接返回
    if boxes:
        return boxes
    
    # 如果沒找到3個數字的框，找4個數字的框
    for r1 in range(rows):
        for c1 in range(cols):
            for r2 in range(r1, rows):
                for c2 in range(c1, cols):
                    # 計算框內格子數
                    cell_count = (r2 - r1 + 1) * (c2 - c1 + 1)
                    # 只考慮4個數字的框
                    if cell_count == 4:
                        if calculate_box_sum(grid, r1, c1, r2, c2) == 10:
                            score = count_nonzero(grid, r1, c1, r2, c2)
                            if score:  # 至少有 1 格還不是 0
                                cells = [(r, c) for r in range(r1, r2 + 1)
                                                 for c in range(c1, c2 + 1)]
                                boxes.append({
                                    'coords': (r1, c1, r2, c2),
                                    'score': score,
                                    'cells': cells
                                })
    
    # 如果找到4個數字的框，直接返回
    if boxes:
        return boxes
    
    # 如果都沒找到，則找所有可能的框
    for r1 in range(rows):
        for c1 in range(cols):
            for r2 in range(r1, rows):
                for c2 in range(c1, cols):
                    # 計算框內格子數
                    cell_count = (r2 - r1 + 1) * (c2 - c1 + 1)
                    # 跳過已經找過的2、3、4個數字的框
                    if cell_count <= 4:
                        continue
                    if calculate_box_sum(grid, r1, c1, r2, c2) == 10:
                        score = count_nonzero(grid, r1, c1, r2, c2)
                        if score:  # 至少有 1 格還不是 0
                            cells = [(r, c) for r in range(r1, r2 + 1)
                                             for c in range(c1, c2 + 1)]
                            boxes.append({
                                'coords': (r1, c1, r2, c2),
                                'score': score,
                                'cells': cells
                            })
    
    return boxes

def solve_greedy(grid):
    """
    使用貪婪演算法來選擇互不重疊的框選，以最大化總分數。
    """
    # 找到所有有效的和為 10 的框選
    valid_boxes = find_valid_boxes(grid)

    # 按分數從高到低排序有效框選
    valid_boxes.sort(key=lambda x: x['score'], reverse=True)

    selected_boxes = [] # 存放最終選擇的框選
    used_cells = set()  # 用於記錄已經被選擇過的單元格座標 (row, col)

    # 遍歷排序後的框選
    for box in valid_boxes:
        is_overlapping = False
        # 檢查當前框選中的任何單元格是否已經被使用
        for r, c in box['cells']:
            if (r, c) in used_cells:
                is_overlapping = True
                break # 只要有一個單元格重疊就退出內層迴圈

        # 如果框選沒有重疊，則選擇它
        if not is_overlapping:
            selected_boxes.append(box)
            # 將這個框選中的所有單元格標記為已使用
            for r, c in box['cells']:
                used_cells.add((r, c))

    # 計算最終的總分數
    total_score = sum(box['score'] for box in selected_boxes)
    return selected_boxes, total_score

def zero_out_box(grid, cells):
    """把框選內的數字全設為 0。"""
    for r, c in cells:
        grid[r][c] = 0

def greedy_with_zeroing(grid):
    """
    一次找一個最佳框選 → 把它歸零 → 再找下一個，
    直到沒有任何『和 = 10』的矩形。
    回傳 (selected_boxes, total_score)。
    """
    grid = [row[:] for row in grid]          # 深拷貝，別改到原陣列
    selected_boxes = []
    total_score = 0

    while True:
        # 1. 找現在網格裡所有有效框選
        valid_boxes = find_valid_boxes(grid)
        if not valid_boxes:                  # 沒得選就結束
            break

        # 2. 挑分數最高（面積最大）的那一個
        best = max(valid_boxes, key=lambda b: b['score'])
        
        # 3. 記錄被消除的非零數字
        eliminated_nums = []
        for r, c in best['cells']:
            if grid[r][c] != 0:
                eliminated_nums.append(grid[r][c])
        
        # 4. 將被消除的非零數字添加到框選信息中
        best['eliminated_numbers'] = eliminated_nums
        
        selected_boxes.append(best)
        total_score += best['score']

        # 5. 把這個框選歸零
        zero_out_box(grid, best['cells'])

    return selected_boxes, total_score

def random_with_zeroing(grid):
    """
    使用隨機策略來選擇框選，每次隨機選擇一個有效框選，然後把它歸零，
    直到沒有任何『和 = 10』的矩形。
    回傳 (selected_boxes, total_score)。
    """
    grid = [row[:] for row in grid]          # 深拷貝，別改到原陣列
    selected_boxes = []
    total_score = 0

    while True:
        # 1. 找現在網格裡所有有效框選
        valid_boxes = find_valid_boxes(grid)
        if not valid_boxes:                  # 沒得選就結束
            break

        # 2. 隨機選擇一個框選
        best = random.choice(valid_boxes)
        
        # 3. 記錄被消除的非零數字
        eliminated_nums = []
        for r, c in best['cells']:
            if grid[r][c] != 0:
                eliminated_nums.append(grid[r][c])
        
        # 4. 將被消除的非零數字添加到框選信息中
        best['eliminated_numbers'] = eliminated_nums
        
        selected_boxes.append(best)
        total_score += best['score']

        # 5. 把這個框選歸零
        zero_out_box(grid, best['cells'])

    return selected_boxes, total_score

def run_multiple_random(grid, num_runs=200, max_workers=20):
    """
    並行執行多次隨機算法，找出得分最高的那次結果。
    
    參數:
    - grid: 初始網格數據
    - num_runs: 執行次數，默認為200
    - max_workers: 最大並行工作數，默認為20
    
    返回:
    - 得分最高的框選列表和總分
    """
    print(f"開始並行執行 {num_runs} 次隨機算法...")
    start_time = time.time()
    
    # 使用進程池並行執行
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        # 提交所有任務
        futures = {executor.submit(random_with_zeroing, grid): i for i in range(num_runs)}
        
        # 收集結果，使用 tqdm 顯示進度
        results = []
        for future in tqdm(as_completed(futures), total=num_runs, desc="執行進度"):
            run_id = futures[future]
            try:
                boxes, score = future.result()
                results.append((boxes, score, run_id))
            except Exception as e:
                print(f"運行 {run_id} 時發生錯誤: {e}")
    
    # 找出得分最高的結果
    best_result = max(results, key=lambda x: x[1])
    # for i in range(len(results)):
    #     print(results[i][1])
    best_boxes, best_score, best_run_id = best_result
    
    end_time = time.time()
    total_time = end_time - start_time
    
    print(f"\n完成 {num_runs} 次隨機算法執行，總耗時: {total_time:.2f} 秒")
    print(f"最高分數: {best_score} (來自第 {best_run_id+1} 次運行)")
    
    return best_boxes, best_score, best_run_id

def visualize_elimination_steps(grid, boxes):
    """
    顯示每次消除後的矩陣，並輸出座標和消除的數字。
    
    參數:
    - grid: 初始網格數據
    - boxes: 消除的框選列表，每個元素包含座標、分數和消除的數字
    """
    # 深拷貝初始網格，避免修改原始數據
    current_grid = [row[:] for row in grid]
    
    print("\n初始矩陣:")
    print_grid(current_grid)
    
    # 遍歷每個框選，顯示消除過程
    for i, box in enumerate(boxes):
        r1, c1, r2, c2 = box['coords']
        eliminated_nums = box['eliminated_numbers']
        score = box['score']
        
        print(f"\n第 {i+1} 步消除:")
        print(f"  座標: ({r1+1}, {c1+1}) 到 ({r2+1}, {c2+1})")
        print(f"  被消除的非零數字: {eliminated_nums}")
        print(f"  分數: {score}")
        
        # 先將所有X標記清除為空白
        for r in range(len(current_grid)):
            for c in range(len(current_grid[0])):
                if current_grid[r][c] == "X":
                    current_grid[r][c] = " "
        
        # 將當前框選內的數字設為X標記
        for r, c in box['cells']:
            current_grid[r][c] = "X"
        
        # 顯示消除後的矩陣
        print("消除後的矩陣:")
        print_grid(current_grid)
    
    print(f"\n總共消除了 {len(boxes)} 個框選，總分數: {sum(box['score'] for box in boxes)}")

def save_elimination_steps_to_file(grid, boxes, filename="elimination_steps.txt"):
    """
    將消除步驟保存到文件中。
    
    參數:
    - grid: 初始網格數據
    - boxes: 消除的框選列表
    - filename: 輸出文件名
    """
    # 深拷貝初始網格，避免修改原始數據
    current_grid = [row[:] for row in grid]
    
    with open(filename, 'w', encoding='utf-8') as f:
        f.write("初始矩陣:\n")
        for row in current_grid:
            f.write("  ".join(str(cell) for cell in row) + "\n")
        
        # 遍歷每個框選，記錄消除過程
        for i, box in enumerate(boxes):
            r1, c1, r2, c2 = box['coords']
            eliminated_nums = box['eliminated_numbers']
            score = box['score']
            
            f.write(f"\n第 {i+1} 步消除:\n")
            f.write(f"  座標: ({r1+1}, {c1+1}) 到 ({r2+1}, {c2+1})\n")
            f.write(f"  被消除的非零數字: {eliminated_nums}\n")
            f.write(f"  分數: {score}\n")
            
            # 先將所有X標記清除為空白
            for r in range(len(current_grid)):
                for c in range(len(current_grid[0])):
                    if current_grid[r][c] == "X":
                        current_grid[r][c] = " "
            
            # 將當前框選內的數字設為X標記
            for r, c in box['cells']:
                current_grid[r][c] = "X"
            
            # 記錄消除後的矩陣
            f.write("消除後的矩陣:\n")
            for row in current_grid:
                for cell in row:
                    if cell == "X":
                        f.write(" - ")
                    else:
                        f.write("  "+str(cell))
                f.write("\n")
                # if row == "X":
                #     f.write(" -".join(str(cell) for cell in row) + "\n")
                # else:
                #     f.write("  ".join(str(cell) for cell in row) + "\n")
        
        f.write(f"\n總共消除了 {len(boxes)} 個框選，總分數: {sum(box['score'] for box in boxes)}\n")
    
    print(f"消除步驟已保存到文件: {filename}")

def print_grid(grid):
    """格式化輸出矩陣"""
    for row in grid:
        print("  ".join(str(cell) for cell in row))

# 根據圖片建立範例網格數據 (10行 x 18列)
# 這裡假設圖片中的空白蘋果位置數字為 0
# grid_data = [
#     [7, 2, 6, 3, 6, 6, 4, 7, 4, 7, 1, 5, 3, 2, 1, 8, 3],
#     [7, 6, 3, 1, 2, 1, 2, 2, 6, 4, 8, 1, 2, 4, 8, 1, 7],
#     [7, 3, 2, 9, 7, 5, 1, 1, 2, 3, 6, 9, 1, 5, 5, 9, 2],
#     [7, 3, 5, 5, 4, 4, 9, 5, 7, 5, 5, 8, 4, 9, 2, 1, 4],
#     [5, 7, 9, 9, 8, 5, 6, 6, 9, 1, 6, 7, 2, 5, 1, 8, 5],
#     [6, 5, 2, 2, 2, 4, 3, 6, 6, 9, 1, 7, 1, 8, 6, 8, 2],
#     [7, 3, 5, 3, 2, 4, 5, 3, 2, 6, 6, 6, 7, 1, 6, 6, 6],
#     [1, 6, 7, 1, 8, 7, 9, 3, 1, 4, 6, 1, 3, 7, 8, 1, 9],
#     [8, 7, 7, 3, 7, 1, 6, 9, 6, 4, 2, 8, 9, 8, 4, 9, 4],
#     [1, 9, 7, 8, 5, 5, 3, 5, 4, 3, 2, 8, 3, 7, 9, 3, 5]
# ]
# grid_data = [
#     [5, 6, 3, 6, 1, 9, 5, 5, 8, 1, 3, 2, 8, 3, 4, 1, 8],
#     [1, 6, 2, 4, 4, 1, 4, 9, 6, 3, 2, 4, 5, 7, 5, 9, 5],
#     [1, 5, 9, 7, 6, 1, 2, 2, 5, 7, 1, 7, 4, 8, 6, 6, 8],
#     [3, 1, 7, 3, 8, 3, 5, 3, 8, 9, 8, 5, 4, 3, 5, 7, 1],
#     [6, 4, 2, 1, 9, 4, 3, 3, 8, 2, 1, 8, 7, 5, 2, 1, 8],
#     [6, 4, 2, 4, 9, 7, 1, 8, 7, 9, 9, 1, 2, 7, 2, 4, 1],
#     [5, 7, 8, 4, 8, 2, 9, 5, 9, 4, 5, 3, 8, 2, 5, 6, 9],
#     [5, 5, 6, 8, 6, 7, 7, 4, 2, 9, 0, 4, 7, 1, 2, 6, 3],
#     [1, 8, 7, 4, 9, 5, 2, 4, 1, 5, 5, 8, 6, 8, 4, 1, 2],
#     [6, 5, 6, 2, 5, 6, 3, 7, 4, 8, 1, 9, 4, 8, 8, 4, 1]
# ]
grid_data = [
    [1, 8, 4    , 4, 2, 1   , 2, 9, 1   , 2, 9, 8   , 1, 2, 7   , 6, 4],
    [4, 2, 2    , 7, 9, 7   , 9, 9, 7   , 5, 4, 5   , 5, 8, 1   , 4, 1],
    [1, 4, 3    , 7, 3, 3   , 5, 6, 3   , 5, 8, 7   , 1, 4, 8   , 9, 8],

    [4, 7, 5    , 9, 2, 6   , 6, 3, 1   , 6, 7, 8   , 1, 6, 9   , 2, 6],
    [8, 7, 3    , 5, 1, 5   , 8, 2, 4   , 1, 5, 1   , 6, 7, 4   , 4, 3],
    [7, 5, 2    , 9, 2, 1   , 6, 7, 9   , 3, 4, 3   , 8, 1, 2   , 4, 3],

    [5, 5, 8    , 5, 1, 1   , 8, 4, 1   , 8, 2, 2   , 2, 3, 1   , 1, 1],
    [7, 5, 2    , 5, 2, 1   , 5, 3, 1   , 9, 1, 2   , 7, 5, 8   , 5, 3],
    [4, 1, 1    , 3, 5, 9   , 2, 9, 6   , 1, 3, 3   , 7, 6, 8   , 9, 9],

    [4, 6, 7    , 4, 4, 9   , 8, 4, 9   , 7, 6, 9   , 3, 5, 3   , 1, 3]
]
grid_data = [
    [7, 4, 9, 2, 1, 6, 7, 1, 5, 2, 9, 1, 4, 3, 8, 6, 1],
    [2, 6, 1, 1, 6, 6, 3, 2, 5, 6, 4, 9, 8, 5, 6, 7, 1],
    [5, 2, 4, 9, 1, 2, 7, 2, 3, 2, 9, 9, 4, 6, 5, 7, 6],
    [4, 5, 4, 9, 6, 9, 8, 3, 6, 1, 3, 2, 1, 9, 2, 3, 5],
    [8, 5, 3, 5, 5, 8, 6, 3, 9, 7, 6, 4, 5, 1, 5, 8, 8],
    [3, 8, 4, 3, 1, 3, 9, 5, 1, 4, 1, 4, 7, 1, 7, 3, 1],
    [8, 9, 1, 9, 7, 6, 3, 8, 4, 7, 5, 9, 7, 4, 9, 4, 2],
    [1, 5, 2, 6, 3, 8, 3, 7, 9, 4, 3, 4, 1, 7, 9, 2, 1],
    [7, 9, 9, 2, 4, 4, 9, 3, 4, 1, 0, 7, 1, 3, 6, 5, 5],
    [3, 2, 2, 4, 6, 1, 2, 7, 2, 9, 8, 5, 4, 8, 4, 3, 5]
]
# 執行貪婪演算法
print("執行貪婪演算法...")
start_time = time.time()
greedy_boxes, greedy_score = greedy_with_zeroing(grid_data)
end_time = time.time()
greedy_time = end_time - start_time

# 執行多次隨機演算法並找出最佳結果
best_boxes, best_score, best_run_id = run_multiple_random(grid_data, num_runs=100000, max_workers=60)

# 輸出貪婪演算法結果
print("\n貪婪演算法選擇的框選：")
# if greedy_boxes:
#     for box in greedy_boxes:
#         print(f"  座標: ({box['coords'][0]+1}, {box['coords'][1]+1}) 到 ({box['coords'][2]+1}, {box['coords'][3]+1}), 分數: {box['score']}")
#         print(f"  被消除的非零數字: {box['eliminated_numbers']}")
# else:
#     print("  沒有找到任何有效且互不重疊的框選。")

print(f"\n貪婪演算法得到的總分數: {greedy_score}")
print(f"貪婪演算法執行時間: {greedy_time:.4f} 秒")

# 輸出最佳隨機演算法結果
print(f"\n最佳隨機演算法結果 (第 {best_run_id+1} 次運行):")
if best_boxes:
    for box in best_boxes:
        # print(f"  座標: ({box['coords'][0]+1}, {box['coords'][1]+1}) 到 ({box['coords'][2]+1}, {box['coords'][3]+1}), 分數: {box['score']}")
        if box['coords'][0] == box['coords'][2]:
            print(f"  座標: ({box['coords'][0]+1}, {box['coords'][1]+1})")
            print(f"  座標: ( , {box['coords'][3]+1})")
        elif box['coords'][1] == box['coords'][3]:
            print(f"  座標: ({box['coords'][0]+1}, {box['coords'][1]+1})")
            print(f"  座標: ({box['coords'][2]+1}, )")
        else:
            print(f"  座標: ({box['coords'][0]+1}, {box['coords'][1]+1})")
            print(f"  座標: ({box['coords'][2]+1}, {box['coords'][3]+1})")
        print(f"  被消除的非零數字: {box['eliminated_numbers']}")
        print(f"  分數: {box['score']}")
else:
    print("  沒有找到任何有效且互不重疊的框選。")

print(f"\n最佳隨機演算法得到的總分數: {best_score}")

# 比較兩種算法的結果
print(f"\n分數差異: {greedy_score - best_score}")
# print(f"貪婪算法 {'優於' if greedy_score > best_score else '劣於'} 最佳隨機算法")

# 顯示最佳隨機算法的消除步驟
print("\n顯示最佳隨機算法的消除步驟:")
# visualize_elimination_steps(grid_data, best_boxes)

# 將消除步驟保存到文件
save_elimination_steps_to_file(grid_data, best_boxes, "/bip8_disk/zhenyu112/bioinfo_handler/handler/fruit_box_elimination_steps.txt")