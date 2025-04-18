#!/usr/bin/env python3
# Example script demonstrating how to use the config parser

import sys
import os

# Add the parent directory to the path so we can import the handler module
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from handler.config_parser import (
    get_config,
    get_config_value,
    get_config_column,
    update_config
)

def main():
    # 基本用法：獲取整個配置
    print("===== 獲取整個配置 =====")
    config = get_config()
    print(f"配置類型: {type(config)}")
    print(f"配置頂層鍵: {list(config.keys())}")
    
    # 獲取特定列
    print("\n===== 獲取特定列 =====")
    data_root_path = get_config_column("data_root_path")
    print(f"data_root_path: {data_root_path}")
    
    # 使用點符號獲取嵌套值
    print("\n===== 使用點符號獲取嵌套值 =====")
    nygc_rank = get_config_value("ans_rank.NYGC")
    print(f"ans_rank.NYGC: {nygc_rank}")
    
    # 使用列表獲取嵌套值
    print("\n===== 使用列表獲取嵌套值 =====")
    seqc2_rank = get_config_value(["ans_rank", "SEQC2"])
    print(f"ans_rank.SEQC2: {seqc2_rank}")
    
    # 使用默認值
    print("\n===== 使用默認值 =====")
    non_existent = get_config_value("non_existent.key", default="默認值")
    print(f"non_existent.key (使用默認值): {non_existent}")
    
    # 獲取列表和字典
    print("\n===== 獲取列表和字典 =====")
    run_samples = get_config_value("run_samples")
    print(f"運行樣本: {run_samples}")
    
    ans_rank = get_config_value("ans_rank")
    print(f"答案排名: {ans_rank}")
    
    # 嘗試更新配置（注意：這會修改實際的配置文件）
    print("\n===== 更新配置（已注釋掉以防止意外修改） =====")
    # update_config({"test_key": "test_value"})
    # print("配置已更新，添加了test_key")
    
    print("\n配置解析器演示完成！")

if __name__ == "__main__":
    main() 