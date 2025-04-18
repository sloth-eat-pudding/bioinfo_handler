from handler.utils import *
import numpy as np
import pandas as pd
import os
import re
import logging
import glob
import joblib
import matplotlib.pyplot as plt
import plotly.graph_objects as go
from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression, Ridge
from sklearn.pipeline import Pipeline
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.model_selection import train_test_split

def extract_features_from_csv(file_path):
    """
    Extract features from CSV file
    
    Args:
        file_path (str): Path to the CSV file
        
    Returns:
        tuple: (features_dict, purity_value, sample_name) or (None, None, None) if error
    """
    try:
        # Read CSV file
        try:
            if file_path.endswith(".q1_q3"):
                # Read file without headers; the first value represents 25% and the second represents 75%
                values = pd.read_csv(file_path, header=None, sep="\t").iloc[0].tolist()
                df = pd.DataFrame(values, index=["25%", "75%"], columns=["ratio"])
            else:
                df = pd.read_csv(file_path, index_col=0)
            # print(df)
        except pd.errors.EmptyDataError:
            logging.error(f"文件 {file_path} 為空，無法處理")
            return None, None, None
        except pd.errors.ParserError:
            logging.error(f"文件 {file_path} 格式錯誤")
            return None, None, None
        
        # Extract features
        try:
            features = {
                # 'mean': df.loc['mean', 'ratio'],
                # 'std': df.loc['std', 'ratio'],
                '25%': df.loc['25%', 'ratio'],
                # '50%': df.loc['50%', 'ratio'],
                '75%': df.loc['75%', 'ratio']
            }
        except KeyError:
            logging.error(f"{file_path} 缺少必要的 `ratio` 統計數據")
            return None, None, None
        
        # Extract purity from filename
        try:
            file_name = os.path.basename(file_path)
            match = re.findall(r't(\d+)_n(\d+)', file_name)
            
            if match:
                t, n = match[0]
                purity = round(int(t) / (int(t) + int(n)), 1)
            else:
                # Try to extract from directory name
                dir_name = os.path.basename(os.path.dirname(file_path))
                match = re.findall(r't(\d+)_n(\d+)', dir_name)
                
                if match:
                    t, n = match[0]
                    purity = round(int(t) / (int(t) + int(n)), 1)
                else:
                    logging.error(f"無法從 {file_path} 提取純度信息")
                    return None, None, None
            
            sample = os.path.basename(os.path.dirname(file_path))
            return features, purity, sample
        except Exception as e:
            logging.error(f"提取純度時出錯: {e}")
            return None, None, None
            
    except Exception as e:
        logging.error(f"處理文件 {file_path} 時出錯: {e}")
        return None, None, None

def collect_csv_files(input_paths, name_pattern="*_1.csv"):
    """
    Collect all CSV files from input paths
    
    Args:
        input_paths (list): List of input file or directory paths
        
    Returns:
        list: List of CSV file paths
    """
    files = []
    for input_path in input_paths:
        if os.path.isfile(input_path):
            files.append(input_path)
        elif os.path.isdir(input_path):
            files.extend(glob.glob(os.path.join(input_path, name_pattern)))
        else:
            raise ValueError(f"輸入路徑 {input_path} 不存在或不是有效的文件/目錄")
    
    if not files:
        raise ValueError("在提供的路徑中沒有找到CSV文件")
    
    return files

def extract_data_from_files(files):
    """
    Extract features and purity from all files
    
    Args:
        files (list): List of CSV file paths
        
    Returns:
        tuple: (X_data, y_data, sample_data, key_data)
    """
    data = [extract_features_from_csv(file) for file in files]
    valid_data = [(f, p, s) for f, p, s in data if f is not None and len(f) > 0 and p is not None]
    
    if not valid_data:
        raise ValueError("沒有找到有效的數據進行訓練")
    
    # Separate features, purity values and sample names
    X_data = [list(f.values()) for f, _, _ in valid_data]
    y_data = [p for _, p, _ in valid_data]
    sample_data = [s for _, _, s in valid_data]
    
    # Get feature names from the first valid sample
    key_data = list(valid_data[0][0].keys())
    
    return X_data, y_data, sample_data, key_data

def train_polynomial_model(X, y, degree=2, model_path=None):
    """
    Train polynomial regression model
    
    Args:
        X (numpy.ndarray): Feature matrix
        y (numpy.ndarray): Target values
        degree (int): Polynomial degree
        model_path (str): Path to save the model
        
    Returns:
        tuple: (model, mse, r2)
    """
    # Create polynomial regression model with Ridge regularization to prevent overfitting
    model = Pipeline([
        ('poly', PolynomialFeatures(degree=degree)),
        # ('ridge', Ridge(alpha=0.1))  # Use Ridge regression with L2 regularization
        ('linear', LinearRegression())
    ])
    
    # Train model on all data
    model.fit(X, y)
    y_pred = model.predict(X)
    
    # Calculate metrics
    mse = mean_squared_error(y, y_pred)
    r2 = r2_score(y, y_pred)
    
    logging.info(f"模型訓練完成 - MSE: {mse:.4f}, R²: {r2:.4f}")
    
    # Save model if path is provided
    if model_path:
        joblib.dump(model, model_path)
        logging.info(f"模型已保存到 {model_path}")
        
        # Get coefficients and intercept
        coefficients = model.named_steps['linear'].coef_
        # coefficients = model.named_steps['ridge'].coef_
        intercept = model.named_steps['linear'].intercept_
        feature_names = model.named_steps['poly'].get_feature_names_out()
        terms = [f"{coef:.4f}*{name}" for coef, name in zip(coefficients, feature_names)]
        equation = f"{intercept:.4f} + " + " + ".join(terms)
        equation = equation.replace("+ -", "- ")
        # equation = f"{intercept:.4f} + " + " + ".join(terms).replace("+ -", "- ")
        logging.info(f"模型方程: f(x) = {equation}")
    
    return model, mse, r2

def visualize_model_3d(X, y_data, predictions, sample_data, key_data, model, model_path=None):
    """
    Create 3D visualization of the model
    
    Args:
        X (numpy.ndarray): Feature matrix
        y_data (list): Actual purity values
        predictions (numpy.ndarray): Predicted purity values
        sample_data (list): Sample names
        key_data (list): Feature names
        model (sklearn.pipeline.Pipeline): Trained model
        model_path (str): Path to save the visualization
    """
    # Create color mapping for unique purity values
    unique_purities = sorted(set(y_data))
    num_classes = len(unique_purities)
    cmap = plt.get_cmap("tab10", num_classes)
    purity_to_color = {p: cmap(i)[:3] for i, p in enumerate(unique_purities)}
    
    # Create 3D prediction surface
    x0 = np.linspace(0, 1, 100)
    x1 = np.linspace(0, 1, 100)
    X0, X1 = np.meshgrid(x0, x1)
    grid = np.c_[X0.ravel(), X1.ravel()]
    
    # Predict values on the grid
    y_pred = model.predict(grid)
    y_pred = np.clip(y_pred, 0, 1)
    Y = y_pred.reshape(X0.shape)
    
    # ======== Matplotlib 3D plot ========
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection='3d')
    surf = ax.plot_surface(X0, X1, Y, cmap='viridis', alpha=0.6)
    
    # Store legend handles
    scatter_plots = {}
    
    # Plot data points
    for i in range(len(X)):
        purity = y_data[i]
        color = purity_to_color[purity]
        error = abs(predictions[i] - purity)
        
        if purity not in scatter_plots:
            scatter_plots[purity] = ax.scatter([], [], [], color=color, label=f"purity {purity}")
        
        # Size based on error (larger = more error)
        # point_size = 40 + error * 200
        point_size = 40
        ax.scatter(X[i, 0], X[i, 1], predictions[i], color=color, s=point_size)
        
        # Add text label with sample name and error
        ax.text(X[i, 0], X[i, 1], predictions[i], 
               f"{sample_data[i]}\nerror: {error:.2%}", 
               fontsize=8, color=color)
    
    # Set axis labels and view
    ax.set_xlabel(key_data[0])
    ax.set_ylabel(key_data[1])
    ax.set_zlabel('purity')
    ax.view_init(elev=30, azim=45)
    
    # Add colorbar for surface
    fig.colorbar(surf, ax=ax, shrink=0.5, aspect=5, label='purity')
    
    # Add legend
    ax.legend(handles=list(scatter_plots.values()), loc='best')
    
    # Save figure
    if model_path:
        plt.savefig(model_path.replace('.joblib', '_prediction.png'), dpi=300, bbox_inches='tight')
    
    plt.show()
    
    # ======== Plotly interactive 3D plot ========
    fig = go.Figure()
    
    # Add surface
    fig.add_trace(go.Surface(z=Y, x=X0, y=X1, colorscale='Viridis', opacity=0.6))
    
    # Plot data points with hover information
    for purity in unique_purities:
        indices = [i for i, p in enumerate(y_data) if p == purity]
        color = f"rgb({int(purity_to_color[purity][0]*255)}, {int(purity_to_color[purity][1]*255)}, {int(purity_to_color[purity][2]*255)})"
        
        # Create hover text with detailed information
        hover_texts = [
            f"樣本: {sample_data[i]}<br>" +
            f"真實純度: {y_data[i]:.2f}<br>" +
            f"預測純度: {predictions[i]:.2f}<br>" +
            f"誤差: {abs(predictions[i] - y_data[i]):.2%}<br>" +
            f"特徵 - {key_data[0]}: {X[i, 0]:.4f}<br>" +
            f"特徵 - {key_data[1]}: {X[i, 1]:.4f}"
            for i in indices
        ]
        mode = f'markers'
        
        # Add scatter trace
        fig.add_trace(go.Scatter3d(
            x=[X[i, 0] for i in indices],
            y=[X[i, 1] for i in indices],
            z=[predictions[i] for i in indices],
            mode=mode,
            marker=dict(size=5, color=color),
            text=hover_texts,
            hoverinfo='text',
            name=f"purity {purity}"
        ))
    
    # Set chart layout
    fig.update_layout(
        title='互動式 3D 純度預測模型',
        scene=dict(
            xaxis_title=key_data[0],
            yaxis_title=key_data[1],
            zaxis_title='purity',
            aspectratio=dict(x=1, y=1, z=0.8)
        ),
        legend=dict(
            yanchor="top",
            y=0.99,
            xanchor="left",
            x=0.01
        ),
        margin=dict(l=0, r=0, b=0, t=30)
    )
    
    # Save as HTML if model path is provided
    if model_path:
        html_path = model_path.replace('.joblib', '_interactive.html')
        fig.write_html(html_path)
        logging.info(f"互動式圖表已保存到 {html_path}")
    
    fig.show()

def print_grouped_data(data):
    """
    Print prediction results grouped by actual purity values
    
    Args:
        data (list): List of tuples (actual_purity, predicted_purity, error, sample_name)
    """
    from collections import defaultdict
    grouped_data = defaultdict(list)

    # 依據 y_data[i] 分組
    for y, pred, err, sample in data:
        grouped_data[y].append((pred, err, sample))

    # 輸出結果
    for key, values in grouped_data.items():
        print(f"actual:{key}")
        print(f"  Predicted, Error , Sample")
        for v in values:
            print(f"  {v[0]:.5f}, {v[1]:.5f}, {v[2]}")
        print()  # 空行增加可讀性

def plot_prediction_vs_actual(prediction_data, output_path='purity_prediction_results.png'):
    """
    Plot the prediction vs actual purity values with different colors for each sample.
    Connect points from the same sample with lines.
    
    Args:
        prediction_data: List of tuples (actual_purity, predicted_purity, error, sample_name)
        output_path: Path to save the output image
    """
    import matplotlib.pyplot as plt
    import numpy as np
    
    # Extract data from tuples
    actual_values = [item[0] for item in prediction_data]
    predicted_values = [item[1] for item in prediction_data]
    sample_names = [item[3] for item in prediction_data]
    
    # Get unique sample names for coloring
    unique_samples = list(set(sample_names))
    # colors = plt.cm.tab10(np.linspace(0, 1, len(unique_samples)))
    
    # Create a color map for samples
    # sample_color_map = {sample: color for sample, color in zip(unique_samples, colors)}
    
    # Create the plot
    plt.figure(figsize=(10, 8))
    
    # Plot perfect prediction line
    plt.plot([0, 1], [0, 1], 'k--', label='Perfect prediction')
    
    # Plot each sample with different color and connect points with lines
    for sample in unique_samples:
        # Get indices for this sample
        sample_indices = [i for i, name in enumerate(sample_names) if name == sample]
        
        # Get actual and predicted values for this sample
        sample_actual = [actual_values[i] for i in sample_indices]
        sample_predicted = [predicted_values[i] for i in sample_indices]
        
        # Sort points by actual purity to ensure proper line connection
        sorted_points = sorted(zip(sample_actual, sample_predicted))
        sorted_actual = [p[0] for p in sorted_points]
        sorted_predicted = [p[1] for p in sorted_points]
        
        # Plot points
        plt.scatter(
            sample_actual,
            sample_predicted,
            color=sample_color_map[sample],
            label=sample,
            alpha=0.8
        )
        
        # Connect points with lines
        plt.plot(
            sorted_actual,
            sorted_predicted,
            color=sample_color_map[sample],
            alpha=0.6,
            linestyle='-'
        )
    
    # Add labels and title
    plt.xlabel('Actual Purity')
    plt.ylabel('Predicted Purity')
    plt.title('Predicted vs Actual Tumor Purity')
    plt.grid(True, linestyle='--', alpha=0.7)
    
    # Add legend
    plt.legend(loc='best')
    
    # Set axis limits
    plt.xlim(0, 1.05)
    plt.ylim(0, 1.05)
    
    # Save the figure
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.show()

# 新增面向对象的类
class PurityPredictor:
    """
    肿瘤純度預測器類
    
    用於訓練和預測肿瘤純度的面向對象實現
    """
    
    def __init__(self, degree=2):
        """
        初始化預測器
        
        Args:
            degree (int): 多項式回歸的階數，默認為2
        """
        self.degree = degree
        self.model = None
        self.model_info = None
        self.feature_names = None
        self.sample_data = None
        self.equation = None  # 存储模型方程
        self.coefficients = None  # 存储模型系数
        self.intercept = None  # 存储模型截距
        self.feature_names_poly = None  # 存储多项式特征名称
        self.predict_result = None
        
    def train(self, input_paths, model_path=None, name_pattern="*.q1_q3"):
        """
        訓練模型
        
        Args:
            input_paths (list): 輸入文件或目錄路徑的列表
            model_path (str): 模型保存路徑，默認為None
            name_pattern (str): 文件名匹配模式，默認為"*.q1_q3"
            
        Returns:
            dict: 包含模型信息和評估指標的字典
        """
        # 收集所有CSV文件
        files = collect_csv_files(input_paths, name_pattern=name_pattern)
        
        # 提取特徵和純度
        X_data, y_data, sample_data, key_data = extract_data_from_files(files)
        
        # 保存特徵名稱和樣本數據
        self.feature_names = key_data
        self.sample_data = sample_data
        
        # 轉換為numpy數組
        X = np.array(X_data)
        y = np.array(y_data)
        
        # 訓練模型
        self.model, mse, r2 = train_polynomial_model(X, y, self.degree, model_path)
        
        # 保存模型信息
        self.model_info = {
            'mse': mse,
            'r2': r2,
            'features': key_data
        }
        
        # 提取并保存模型方程
        if self.model is not None:
            self.coefficients = self.model.named_steps['linear'].coef_
            self.intercept = self.model.named_steps['linear'].intercept_
            self.feature_names_poly = self.model.named_steps['poly'].get_feature_names_out()
            # print(self.feature_names_poly)
            
            # 构建方程字符串
            terms = [f"{coef:.4f}*{name}" for coef, name in zip(self.coefficients, self.feature_names_poly)]
            self.equation = f"{self.intercept:.4f} + " + " + ".join(terms)
            self.equation = self.equation.replace("+ -", "- ")
            
            logging.info(f"模型方程: f(x) = {self.equation}")
        
        return self.model_info
    
    def predict(self, input_paths, model_path=None, plot_model=False, use_function=False, name_pattern="*.q1_q3"):
        """
        使用模型預測肿瘤純度
        
        Args:
            input_paths (list): 輸入文件或目錄路徑的列表
            model_path (str): 模型路徑，如果為None則使用已加載的模型
            plot_model (bool): 是否繪製模型可視化圖，默認為False
            use_function (bool): 是否使用硬編碼的函數計算而不是模型預測，默認為False
            name_pattern (str): 文件名匹配模式，默認為"*.q1_q3"
            
        Returns:
            float/list: 預測的肿瘤純度值或結果列表
            dict: 包含模型信息和評估指標的字典
        """
        # 收集所有CSV文件
        files = collect_csv_files(input_paths, name_pattern=name_pattern)
        
        # 提取特徵和純度
        X_data, y_data, sample_data, key_data = extract_data_from_files(files)
        
        # 保存特徵名稱和樣本數據
        self.feature_names = key_data
        self.sample_data = sample_data
        
        # 轉換為numpy數組
        X = np.array(X_data)
        y = np.array(y_data)
        
        # 只啟用3D繪圖，如果我們有恰好2個特徵
        plot_model = plot_model and len(key_data) == 2
        
        # 如果使用硬編碼的函數計算而不是模型預測
        if use_function:
            # 使用硬編碼的函數計算純度
            # 這個函數是從 PhasingGraph.cpp 中的 getPurity 函數移植過來的
            # def predict_with_values(q1, q3):
            #     purity = -36.2584 + 0.0000*1 + 178.0867*q1 - 13.2564*q3 + 130.8341*q1*q1 - 559.6774*q1*q3 + 213.7812*q3*q3 - 89.8341*q1*q1*q1 + 5.4682*q1*q1*q3 + 328.7876*q1*q3*q3 - 162.3958*q3*q3*q3
            #     # 確保純度值在 0 和 1 之間
            #     return max(0.0, min(purity, 1.0))
            
            # 如果只有一個文件，直接計算
            # if len(files) == 1:
            #     q1, q3 = X_data[0]
            #     prediction = self.predict_with_values(q1, q3)
            #     logging.info(f"使用函數計算的肿瘤純度: {prediction:.4f}")
            #     return prediction, {'function': 'predict_with_values', 'prediction': prediction}
            
            # 對於多個文件，返回所有計算結果
            predictions = [self.predict_with_values(q1, q3) for q1, q3 in X_data]
            results = []
            
            # 計算每個預測的誤差
            errors = [abs(y[i] - predictions[i]) for i in range(len(predictions))]
            
            for i in range(len(files)):
                results.append((y_data[i], round(predictions[i], 5), round(errors[i], 5), sample_data[i]))
            
            # 根據實際純度值排序結果
            results.sort()
            
            # return results, {'function': 'calculate_purity_function'}
            return results
        else:
            # 使用模型預測
            # 加載模型
            if model_path:
                if not os.path.exists(model_path):
                    raise ValueError("需要提供有效的模型路徑")
                self.model = joblib.load(model_path)
            elif self.model is None:
                raise ValueError("需要提供有效的模型路徑或先訓練模型")
        
        # 如果只有一個文件，直接預測
        # if len(files) == 1:
        #     prediction = self.model.predict([X_data[0]])[0]
        #     prediction = np.clip(prediction, 0, 1)  # 確保預測值在0和1之間
        #     logging.info(f"預測的肿瘤純度: {prediction:.4f}")
        #     return prediction, {'model': self.model, 'prediction': prediction}
        
        # 對於多個文件，返回所有預測結果
        predictions = np.clip(self.model.predict(X), 0, 1)
        results = []
        
        # 計算每個預測的誤差
        errors = [abs(y[i] - predictions[i]) for i in range(len(predictions))]
        
        for i in range(len(files)):
            results.append((y_data[i], round(predictions[i], 5), round(errors[i], 5), sample_data[i]))
        
        # 根據實際純度值排序結果
        results.sort()
        self.predict_result = results
        
        # # 如果啟用，並且我們有2個特徵，則可視化
        if plot_model:
            visualize_model_3d(X, y_data, predictions, sample_data, key_data, self.model, model_path)
        
        return results
    
    def load_model(self, model_path):
        """
        加載已訓練的模型
        
        Args:
            model_path (str): 模型路徑
            
        Returns:
            bool: 是否成功加載模型
        """
        if not os.path.exists(model_path):
            logging.error(f"模型路徑 {model_path} 不存在")
            return False
        
        try:
            self.model = joblib.load(model_path)
            logging.info(f"成功加載模型: {model_path}")
            return True
        except Exception as e:
            logging.error(f"加載模型時出錯: {e}")
            return False
    
    def save_model(self, model_path):
        """
        保存模型
        
        Args:
            model_path (str): 模型保存路徑
            
        Returns:
            bool: 是否成功保存模型
        """
        if self.model is None:
            logging.error("沒有可保存的模型，請先訓練模型")
            return False
        
        try:
            joblib.dump(self.model, model_path)
            logging.info(f"模型已保存到 {model_path}")
            return True
        except Exception as e:
            logging.error(f"保存模型時出錯: {e}")
            return False
    
    def visualize(self, X, y_data, predictions, sample_data, key_data, model_path=None):
        """
        可視化模型
        
        Args:
            X (numpy.ndarray): 特徵矩陣
            y_data (list): 實際純度值
            predictions (numpy.ndarray): 預測純度值
            sample_data (list): 樣本名稱
            key_data (list): 特徵名稱
            model_path (str): 模型路徑，用於保存可視化結果
        """
        if len(key_data) != 2:
            logging.warning("只有當特徵數量為2時才能進行3D可視化")
            return
        
        visualize_model_3d(X, y_data, predictions, sample_data, key_data, self.model, model_path)
    
    def plot_results(self, results, output_path='purity_prediction_results.png'):
        """
        繪製預測結果
        
        Args:
            results (list): 預測結果列表
            output_path (str): 輸出圖像路徑
        """
        plot_prediction_vs_actual(results, output_path)
    
    def print_results(self, results):
        """
        打印預測結果
        
        Args:
            results (list): 預測結果列表
        """
        print_grouped_data(results)

    def make_feats(self, x0, x1, feature_names):
        """
        根据 feature_names 里每个字符串，计算对应的数值特征，
        并返回一个 numpy 向量，顺序和 feature_names 保持一致。
        """
        # 把 x0, x1 放到 eval 的局部环境里
        env = {'x0': x0, 'x1': x1}
        feats = []
        for name in feature_names:
            if name == '1':
                # 常数项直接 1
                feats.append(1)
            else:
                # 1) '^' -> '**'；2) ' ' -> '*'
                expr = name.replace('^', '**').replace(' ', '*')
                # 在受限环境里 eval 出结果
                feats.append(eval(expr, {}, env))
        return np.array(feats)
    def predict_with_values(self, x0, x1):
        """
        使用两个输入变量直接预测肿瘤纯度
        
        Args:
            x1 (float): 第一个特征值
            x2 (float): 第二个特征值
            
        Returns:
            float: 预测的肿瘤纯度值
        """
        if self.model is None:
            raise ValueError("模型未训练，请先训练模型")
            
        if len(self.feature_names) != 2:
            raise ValueError("此方法仅支持两个特征的模型")
        # 使用模型预测
        # prediction = self.model.predict([[x1, x2]])[0]
        # print(self.feature_names_poly)
        feats = self.make_feats(x0, x1, self.feature_names_poly)
        prediction = self.intercept + sum(c * f for c, f in zip(self.coefficients, feats))
        prediction = np.clip(prediction, 0, 1)  # 确保预测值在0和1之间
        
        # logging.info(f"输入值 ({x1}, {x2}) 的预测肿瘤纯度: {prediction:.4f}")
        return prediction
        
    def get_equation(self):
        """
        获取模型方程
        
        Returns:
            str: 模型方程字符串
        """
        if self.equation is None:
            return "模型未训练，无法获取方程"
        return self.equation

# 保留原有的函數，但標記為已棄用
def predict_purity_polynomial_regression(input_paths, degree=2, train_model=True, model_path=None, plot_model=False, use_function=False):
    """
    使用多項式回歸預測肿瘤純度 (已棄用，請使用 PurityPredictor 類)
    
    參數:
    input_paths (list): 輸入文件或目錄路徑的列表
    degree (int): 多項式回歸的階數，默認為2
    train_model (bool): 是否訓練新模型，默認為True
    model_path (str): 模型保存路徑，默認為None
    plot_model (bool): 是否繪製模型可視化圖，默認為False
    use_function (bool): 是否使用硬編碼的函數計算而不是模型預測，默認為False
    
    返回:
    float/list: 預測的肿瘤純度值或結果列表
    dict: 包含模型信息和評估指標的字典
    """
    import warnings
    warnings.warn("此函數已棄用，請使用 PurityPredictor 類", DeprecationWarning)
    
    predictor = PurityPredictor(degree=degree)
    
    if train_model:
        return None, predictor.train(input_paths, model_path)
    else:
        return predictor.predict(input_paths, model_path, plot_model, use_function)

