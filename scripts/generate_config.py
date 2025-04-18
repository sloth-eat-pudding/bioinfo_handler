#!/usr/bin/env python3
# Command-line script to generate a default config file

import os
import sys
import argparse

# Add the parent directory to the path so we can import the handler module
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from handler.config_parser import generate_default_config

def main():
    parser = argparse.ArgumentParser(description='生成預設配置文件')
    parser.add_argument('-o', '--output', type=str, help='配置文件輸出路徑')
    parser.add_argument('-f', '--force', action='store_true', help='強制覆蓋現有文件')
    args = parser.parse_args()
    
    try:
        if args.output:
            config_path = generate_default_config(output_path=args.output, overwrite=args.force)
        else:
            config_path = generate_default_config(overwrite=args.force)
        
        print(f"成功生成預設配置文件：{config_path}")
        print("您現在可以編輯此文件以適應您的需求。")
        
    except FileExistsError as e:
        print(f"錯誤：{e}")
        print("提示：使用 -f 或 --force 參數可以覆蓋現有文件")
        sys.exit(1)
    except PermissionError as e:
        print(f"錯誤：{e}")
        print("提示：請確保您有寫入權限")
        sys.exit(1)
    except Exception as e:
        print(f"發生未知錯誤：{e}")
        sys.exit(1)

if __name__ == "__main__":
    main() 