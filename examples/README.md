# 示例脚本

本目录包含用于测试和演示`bioinfo_handler`库功能的示例脚本。

## 测试PrecisionCalculator的output_picture方法

`test_output_picture.py`脚本演示了如何使用`PrecisionCalculator`类的`output_picture`方法生成雷达图，用于可视化不同样本、纯度和软件的性能指标（Precision、Recall、F1-score）。

### 使用方法

1. 确保已安装所需依赖：
   ```bash
   pip install pandas matplotlib numpy tqdm
   ```

2. 运行测试脚本：
   ```bash
   python examples/test_output_picture.py
   ```

3. 查看生成的图表：
   - 单独的雷达图位于`output_test/metrics_plot/`目录下
   - 合并的雷达图位于`output_test/`目录下

### 示例数据

示例数据存储在`example_data.py`文件中，包含了不同样本（H2009、HCC1954等）、不同纯度（0.2-1.0）和不同软件（DeepSomatic、Longphase、ClairS）的性能指标。

这些数据可以用于测试和演示`bioinfo_handler`库的功能，也可以作为开发新功能的参考。 