# Bioinfo Handler

生物信息學數據處理工具集。

希望它是你開發的好助手 :)

## 主要功能

- **主要執行 (Main)**: 提供核心功能和工作流程執行：
  - 文件路徑管理 (file_for): 智能定位和管理樣本、純度和軟件的文件路徑
  - 批量處理: 支持多樣本、多純度和多軟件的批量處理
  - 數據分析流程: 自動化執行從數據讀取到結果生成的完整流程
  - 結果可視化: 生成盒形圖等數據可視化結果

- **命令執行器 (Command Runner)**: 提供執行生物信息學分析命令的工具，包括：
  - Phasing: 使用 longphase 進行基因型分型
  - Tagging: 為 BAM 文件添加標籤
  - 索引: 為 BAM 文件創建索引
  - 合併: 合併多個文件的結果
  - 支持 Linux 命令執行和日誌記錄

- **文件路徑管理 (File For)**: 提供快速獲取特定樣本、純度和軟件的文件路徑功能：
  - 自動選擇最佳參考數據
  - 根據樣本和純度找到對應的 BAM 文件
  - 根據軟件類型選擇正確的 VCF 輸出路徑
  - 支持多種變異類型 (SNV, Indel)

- **生物信息讀取器 (Bioinfo Reader)**: 提供讀取和處理各種生物信息學文件格式的功能：
  - VCF 文件讀取和解析
  - BED 文件處理
  - 染色體位置排序
  - 變異類型檢查
  - 格式化數據提取

- **配置解析器 (Config Parser)**: 提供靈活的 YAML 配置文件解析和管理功能，支持嵌套配置、默認值和配置更新。

- **純度預測 (Purity Prediction)**: 提供腫瘤樣本純度預測和分析工具。

- **指標計算 (Metrics)**: 提供計算和評估生物信息學分析結果的各種指標。

## 安裝環境

為了安裝和運行此工具集，您需要配置相應的環境。您可以使用以下 `file.yaml` 文件來設置所需的依賴項：

```yaml
dependencies:
  - python=3.8
  - pyyaml
  - pandas
  - matplotlib
  - numpy
  -其他必要的庫
```

您可以使用以下命令來創建環境：

```bash
conda env create -f file.yaml
```

然後，激活環境：

```bash
conda activate your_environment_name
```

確保您已安裝所有必要的依賴項，以便順利運行此工具集。

## 配置解析器 (Config Parser)

這個模塊提供了一個靈活的YAML配置文件解析器，可以輕鬆讀取和更新配置。

### 開始使用前

在使用此工具之前，您需要先創建一個配置文件。您可以通過以下三種方式之一來創建：

#### 1. 使用命令行工具生成預設配置文件

最簡單的方法是使用內置的命令行工具：

```bash
# 生成預設配置文件到默認位置
python -m bioinfo_handler.scripts.generate_config

# 指定自定義輸出路徑
python -m bioinfo_handler.scripts.generate_config -o /path/to/your/custom_config.yaml

# 強制覆蓋現有文件
python -m bioinfo_handler.scripts.generate_config -f
```

命令行參數：
- `-o, --output`: 指定配置文件的輸出路徑
- `-f, --force`: 強制覆蓋現有文件

#### 2. 使用Python API生成預設配置文件

```python
from bioinfo_handler.handler.config_parser import generate_default_config

# 生成預設配置文件到默認位置
config_path = generate_default_config()
print(f"配置文件已生成：{config_path}")

# 或者指定自定義路徑
custom_path = "/path/to/your/custom_config.yaml"
config_path = generate_default_config(output_path=custom_path)
```

您也可以直接運行示例腳本來生成預設配置：

```bash
python -m bioinfo_handler.examples.generate_config_example
```

#### 3. 手動創建配置文件

您可以手動創建一個名為 `config.yaml` 的文件，並放置在項目根目錄下。配置文件應包含以下基本結構：

```yaml
data_root_path: "/your/data/path"
run_root_path: "/your/output/path"
run_samples:
  - sample1
  - sample2
```

### 基本用法

```python
from bioinfo_handler.handler.config_parser import get_config, get_config_value

# 獲取整個配置
config = get_config()

# 獲取特定配置值（使用點符號）
data_root_path = get_config_value("data_root_path")
nygc_rank = get_config_value("ans_rank.NYGC")

# 使用默認值（當配置項不存在時）
timeout = get_config_value("connection.timeout", default=30)
```

### 主要功能

- **generate_default_config(output_path=None, overwrite=False)**: 生成預設配置文件
- **get_config(config_path=None)**: 讀取並返回整個配置字典
- **get_config_value(key_path, default=None, config_path=None)**: 使用點符號或鍵列表獲取嵌套配置值
- **get_config_column(column_name, config_path=None)**: 獲取頂層配置項
- **update_config(updates, config_path=None)**: 更新配置並寫回文件

### 配置文件格式

配置文件使用YAML格式，支持嵌套結構、列表和基本數據類型：

```yaml
data_root_path: "/big8_disk/data"
run_root_path: "/bip8_disk/user"

# 嵌套字典
ans_rank:
  NYGC: 2
  SEQC2: 2

# 列表
run_samples:
  - H2009
  - HCC1954
  - HCC1937
```

### 示例

- 查看 `examples/config_example.py` 獲取配置讀取示例
- 查看 `examples/generate_config_example.py` 獲取配置生成示例

## 安裝依賴

```bash
pip install pyyaml
``` 