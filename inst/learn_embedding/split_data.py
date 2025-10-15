#!/usr/bin/env python3
"""
Split the all_labeled.csv data into train/test sets and save separate CSV files.
This creates:
- train_labeled.csv: 50% of data for training (will be further split into train/val)
- test_labeled.csv: 50% of data for final test/embedding extraction
"""

import argparse
import pandas as pd
from pathlib import Path
from sklearn.model_selection import train_test_split


def split_data(input_csv, output_dir, test_size=0.5, stratify=False, seed=1337):
    """
    Split the labeled data into train and test sets.
    
    Args:
        input_csv: Path to the all_labeled.csv file
        output_dir: Directory to save the split CSV files
        test_size: Fraction of data to use for test set (default: 0.5)
        stratify: Whether to stratify split (for classification tasks)
        seed: Random seed for reproducibility
    """
    # Read the data
    df = pd.read_csv(input_csv)
    print(f"📊 Loaded data: {len(df)} samples")
    print(f"   Columns: {df.columns.tolist()}")
    
    # Check data distribution
    if 'y' in df.columns:
        print(f"   Target stats: mean={df['y'].mean():.2f}, std={df['y'].std():.2f}")
        print(f"   Target range: [{df['y'].min():.2f}, {df['y'].max():.2f}]")
    
    # Determine stratification
    stratify_col = None
    if stratify and 'y' in df.columns:
        # For classification, we might want to stratify
        # Check if y looks like categorical data (small number of unique values)
        unique_vals = df['y'].nunique()
        if unique_vals <= 20:  # Arbitrary threshold for categorical
            print(f"   Stratifying by y ({unique_vals} unique values)")
            stratify_col = df['y']
        else:
            print(f"   Not stratifying (y has {unique_vals} unique values, likely continuous)")
    
    # Split the data
    train_df, test_df = train_test_split(
        df, 
        test_size=test_size, 
        random_state=seed,
        stratify=stratify_col
    )
    
    print(f"📈 Split results:")
    print(f"   Train: {len(train_df)} samples ({len(train_df)/len(df)*100:.1f}%)")
    print(f"   Test:  {len(test_df)} samples ({len(test_df)/len(df)*100:.1f}%)")
    
    # Create output directory
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Save the split files
    train_csv = output_dir / "train_labeled.csv"
    test_csv = output_dir / "test_labeled.csv"
    
    train_df.to_csv(train_csv, index=False)
    test_df.to_csv(test_csv, index=False)
    
    print(f"💾 Saved files:")
    print(f"   Train: {train_csv}")
    print(f"   Test:  {test_csv}")
    
    return train_csv, test_csv


def main():
    parser = argparse.ArgumentParser(description="Split labeled data into train/test sets")
    parser.add_argument("--input_csv", required=True, 
                       help="Path to all_labeled.csv file")
    parser.add_argument("--output_dir", required=True,
                       help="Directory to save split CSV files")
    parser.add_argument("--test_size", type=float, default=0.5,
                       help="Fraction of data for test set (default: 0.5)")
    parser.add_argument("--stratify", action="store_true",
                       help="Stratify split by target variable")
    parser.add_argument("--seed", type=int, default=1337,
                       help="Random seed for reproducibility")
    
    args = parser.parse_args()
    
    print("🔀 Splitting data into train/test sets...")
    train_csv, test_csv = split_data(
        args.input_csv,
        args.output_dir, 
        args.test_size,
        args.stratify,
        args.seed
    )
    print("✅ Data splitting completed!")


if __name__ == "__main__":
    main()
