#!/usr/bin/env python3
"""
Simple script to analyze PyTorch Lightning training results
Works with the automatically generated CSV files
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse
import os

def load_training_data(csv_path):
    """Load and clean the Lightning CSV data"""
    df = pd.read_csv(csv_path)
    
    # Extract epoch-level metrics (non-NaN values)
    epoch_data = df[df['train_loss_epoch'].notna() | df['val_loss_epoch'].notna()].copy()
    
    # Get unique epochs and their metrics
    epochs = []
    train_losses = []
    val_losses = []
    loss_diffs = []
    
    for epoch in epoch_data['epoch'].unique():
        if pd.isna(epoch):
            continue
        epoch_rows = epoch_data[epoch_data['epoch'] == epoch]
        
        # Get train loss for this epoch
        train_loss_row = epoch_rows[epoch_rows['train_loss_epoch'].notna()]
        val_loss_row = epoch_rows[epoch_rows['val_loss_epoch'].notna()]
        loss_diff_row = epoch_rows[epoch_rows['loss_diff'].notna()]
        
        if not train_loss_row.empty and not val_loss_row.empty:
            epochs.append(int(epoch))
            train_losses.append(train_loss_row['train_loss_epoch'].iloc[0])
            val_losses.append(val_loss_row['val_loss_epoch'].iloc[0])
            
            if not loss_diff_row.empty:
                loss_diffs.append(loss_diff_row['loss_diff'].iloc[0])
            else:
                loss_diffs.append(val_losses[-1] - train_losses[-1])
    
    return pd.DataFrame({
        'epoch': epochs,
        'train_loss': train_losses,
        'val_loss': val_losses,
        'loss_diff': loss_diffs
    })

def analyze_training(df):
    """Print training analysis"""
    print("="*60)
    print("TRAINING ANALYSIS")
    print("="*60)
    
    print(f"Total epochs: {len(df)}")
    print(f"Best training loss: {df['train_loss'].min():.4f} (epoch {df.loc[df['train_loss'].idxmin(), 'epoch']})")
    print(f"Best validation loss: {df['val_loss'].min():.4f} (epoch {df.loc[df['val_loss'].idxmin(), 'epoch']})")
    print(f"Final training loss: {df['train_loss'].iloc[-1]:.4f}")
    print(f"Final validation loss: {df['val_loss'].iloc[-1]:.4f}")
    
    # Overfitting analysis
    final_gap = df['loss_diff'].iloc[-1]
    print(f"\nOverfitting Analysis:")
    print(f"Final validation-training gap: {final_gap:.4f}")
    if final_gap > df['train_loss'].iloc[-1]:
        print("⚠️  Significant overfitting detected!")
    elif final_gap > 0.1 * df['train_loss'].iloc[-1]:
        print("⚠️  Mild overfitting detected")
    else:
        print("✅ No significant overfitting")
    
    # Training progress
    if len(df) > 1:
        train_improvement = (df['train_loss'].iloc[0] - df['train_loss'].iloc[-1]) / df['train_loss'].iloc[0] * 100
        val_improvement = (df['val_loss'].iloc[0] - df['val_loss'].iloc[-1]) / df['val_loss'].iloc[0] * 100
        print(f"\nTraining Improvement:")
        print(f"Training loss improved by: {train_improvement:.1f}%")
        print(f"Validation loss improved by: {val_improvement:.1f}%")

def plot_training_curves(df, save_path=None):
    """Create training curves plot"""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    
    # Loss curves
    ax1.plot(df['epoch'], df['train_loss'], 'b-o', label='Training Loss', linewidth=2, markersize=6)
    ax1.plot(df['epoch'], df['val_loss'], 'r-s', label='Validation Loss', linewidth=2, markersize=6)
    ax1.set_xlabel('Epoch')
    ax1.set_ylabel('Loss')
    ax1.set_title('Training & Validation Loss')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Overfitting indicator
    ax2.plot(df['epoch'], df['loss_diff'], 'g-^', linewidth=2, markersize=6)
    ax2.axhline(y=0, color='k', linestyle='--', alpha=0.5)
    ax2.set_xlabel('Epoch')
    ax2.set_ylabel('Validation Loss - Training Loss')
    ax2.set_title('Overfitting Indicator')
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f"Plot saved to: {save_path}")
    
    return fig

def main():
    parser = argparse.ArgumentParser("Analyze PyTorch Lightning training results")
    parser.add_argument("csv_path", help="Path to metrics.csv file")
    parser.add_argument("--save_plot", help="Path to save plot (optional)")
    parser.add_argument("--show_plot", action="store_true", help="Show plot interactively")
    
    args = parser.parse_args()
    
    if not os.path.exists(args.csv_path):
        print(f"Error: CSV file not found: {args.csv_path}")
        return
    
    # Load and analyze data
    print(f"Loading data from: {args.csv_path}")
    df = load_training_data(args.csv_path)
    
    print(f"\nEpoch-level data:")
    print(df)
    
    # Analysis
    analyze_training(df)
    
    # Plot
    fig = plot_training_curves(df, args.save_plot)
    
    if args.show_plot:
        plt.show()
    
    print("\n✅ Analysis complete!")

if __name__ == "__main__":
    main()
