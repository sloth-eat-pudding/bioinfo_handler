#!/usr/bin/env python3
# Example script demonstrating how to generate a default config file

import sys
import os

# Add the parent directory to the path so we can import the handler module
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from handler.config_parser import generate_default_config

def main():
    # 生成預設配置文件
    try:
        # 使用默認路徑生成配置文件
        config_path = generate_default_config()
        print(f"成功生成預設配置文件：{config_path}")
        
        # 如果要指定路徑，可以這樣做：
        # custom_path = "/path/to/your/custom_config.yaml"
        # config_path = generate_default_config(output_path=custom_path)
        # print(f"成功生成預設配置文件到自定義路徑：{custom_path}")
        
        # 如果要覆蓋現有文件，可以這樣做：
        # config_path = generate_default_config(overwrite=True)
        # print(f"成功覆蓋現有配置文件：{config_path}")
        
    except FileExistsError as e:
        print(f"錯誤：{e}")
        print("提示：使用 overwrite=True 參數可以覆蓋現有文件")
    except PermissionError as e:
        print(f"錯誤：{e}")
        print("提示：請確保您有寫入權限")
    
    print("\n預設配置文件生成示例完成！")

if __name__ == "__main__":
    main() 