#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Test script for PrecisionCalculator's output_picture method
This script uses the example data to test the visualization output
"""

import os
import sys
import pandas as pd
from io import StringIO
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm
import logging

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Add parent directory to path to import bioinfo_handler
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from handler.metrics import PrecisionCalculator
from examples.example_data import get_example_data

def main():
    # Create output directory if it doesn't exist
    output_dir = "output_test"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Initialize PrecisionCalculator
    calculator = PrecisionCalculator()
    
    # Get example data
    data = get_example_data()
    
    # Load data into DataFrame
    df = pd.read_csv(StringIO(data), sep="\t")
    
    # Convert to desired format and load into calculator
    calculator.results = pd.DataFrame([{
        "Category": row["Category"],
        "Precision": row["Precision"],
        "Recall": row["Recall"],
        "F1-score": row["F1-score"],
        "TP": row["TP"],
        "FP": row["FP"],
        "FN": row["FN"]
    } for _, row in df.iterrows()])
    
    # Split Category column into sample, purity, software, and variant_type
    format_columns = ["sample", "purity", "software", "variant_type"]
    calculator.results[format_columns] = calculator.results["Category"].str.split(",", expand=True)
    
    # Convert numeric columns to float
    for col in ["Precision", "Recall", "F1-score"]:
        calculator.results[col] = calculator.results[col].astype(float)
    
    # Convert count columns to int
    for col in ["TP", "FP", "FN"]:
        calculator.results[col] = calculator.results[col].astype(int)
    
    # Print the first few rows to verify
    print("Data loaded successfully:")
    print(calculator.results.head())
    
    # Call output_picture method
    print("Generating radar charts...")
    calculator.output_radar_picture(output_dir)
    print(f"Charts have been saved to {output_dir}/metrics_plot/ and {output_dir}/")

if __name__ == "__main__":
    main() 