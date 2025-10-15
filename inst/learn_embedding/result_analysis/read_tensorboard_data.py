#!/usr/bin/env python3
"""
Script to read and analyze TensorBoard data programmatically
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from tensorboard.backend.event_processing.event_accumulator import EventAccumulator
import argparse

def read_tensorboard_logs(log_dir):
    """Read TensorBoard event files and extract metrics"""
    
    # Find the latest version directory
    versions = [d for d in os.listdir(log_dir) if d.startswith('version_')]
    if not versions:
        raise ValueError(f"No version directories found in {log_dir}")
    
    latest_version = sorted(versions, key=lambda x: int(x.split('_')[1]))[-1]
    event_dir = os.path.join(log_dir, latest_version)
    
    print(f"Reading from: {event_dir}")
    
    # Load the event accumulator
    event_acc = EventAccumulator(event_dir)
    event_acc.Reload()
    
    # Get available tags
    scalar_tags = event_acc.Tags()['scalars']
    print(f"Available metrics: {scalar_tags}")
    
    # Extract data
    data = {}
    for tag in scalar_tags:
        scalar_events = event_acc.Scalars(tag)
        steps = [event.step for event in scalar_events]
        values = [event.value for event in scalar_events]
        data[tag] = {'steps': steps, 'values': values}
    
    return data

def plot_training_curves(data, save_path=None):
    """Create training curves from TensorBoard data"""
    
    fig, axes = plt.subplots(2, 2, figsize=(15, 10))
    fig.suptitle('Training Progress from TensorBoard Data', fontsize=16)
    
    # Training and validation loss
    if 'train_loss_epoch' in data and 'val_loss_epoch' in data:
        ax = axes[0, 0]
        ax.plot(data['train_loss_epoch']['steps'], data['train_loss_epoch']['values'], 
                'b-', label='Train Loss', linewidth=2)
        ax.plot(data['val_loss_epoch']['steps'], data['val_loss_epoch']['values'], 
                'r-', label='Val Loss', linewidth=2)
        ax.set_xlabel('Epoch')
        ax.set_ylabel('Loss')
        ax.set_title('Training & Validation Loss')
        ax.legend()
        ax.grid(True, alpha=0.3)
    
    # Loss difference (overfitting indicator)
    if 'loss_diff' in data:
        ax = axes[0, 1]
        ax.plot(data['loss_diff']['steps'], data['loss_diff']['values'], 
                'g-', linewidth=2)
        ax.axhline(y=0, color='k', linestyle='--', alpha=0.5)
        ax.set_xlabel('Epoch')
        ax.set_ylabel('Val Loss - Train Loss')
        ax.set_title('Overfitting Indicator')
        ax.grid(True, alpha=0.3)
    
    # Step-level training loss
    if 'train_loss_step' in data:
        ax = axes[1, 0]
        ax.plot(data['train_loss_step']['steps'], data['train_loss_step']['values'], 
                'b-', alpha=0.7, linewidth=1)
        ax.set_xlabel('Training Step')
        ax.set_ylabel('Loss')
        ax.set_title('Training Loss (Step Level)')
        ax.grid(True, alpha=0.3)
    
    # Step-level validation loss
    if 'val_loss_step' in data:
        ax = axes[1, 1]
        ax.plot(data['val_loss_step']['steps'], data['val_loss_step']['values'], 
                'r-', alpha=0.7, linewidth=1)
        ax.set_xlabel('Validation Step')
        ax.set_ylabel('Loss')
        ax.set_title('Validation Loss (Step Level)')
        ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f"Plot saved to: {save_path}")
    
    return fig

def export_to_csv(data, save_path):
    """Export TensorBoard data to CSV format"""
    
    # Create a comprehensive DataFrame
    all_data = []
    
    for metric_name, metric_data in data.items():
        for step, value in zip(metric_data['steps'], metric_data['values']):
            all_data.append({
                'metric': metric_name,
                'step': step,
                'value': value
            })
    
    df = pd.DataFrame(all_data)
    df.to_csv(save_path, index=False)
    print(f"Data exported to CSV: {save_path}")
    
    # Also create a wide format for easier analysis
    wide_df = df.pivot(index='step', columns='metric', values='value')
    wide_save_path = save_path.replace('.csv', '_wide.csv')
    wide_df.to_csv(wide_save_path)
    print(f"Wide format exported to: {wide_save_path}")
    
    return df, wide_df

def print_summary_stats(data):
    """Print summary statistics"""
    print("\n" + "="*50)
    print("TRAINING SUMMARY STATISTICS")
    print("="*50)
    
    for metric_name, metric_data in data.items():
        values = metric_data['values']
        if values:
            print(f"\n{metric_name}:")
            print(f"  Final value: {values[-1]:.4f}")
            print(f"  Best value:  {min(values):.4f}")
            print(f"  Mean:        {np.mean(values):.4f}")
            print(f"  Std:         {np.std(values):.4f}")

def main():
    parser = argparse.ArgumentParser("Read TensorBoard logs")
    parser.add_argument("--log_dir", required=True, help="TensorBoard logs directory")
    parser.add_argument("--output_dir", default=".", help="Output directory for plots/CSV")
    parser.add_argument("--show_plot", action="store_true", help="Show plot interactively")
    
    args = parser.parse_args()
    
    # Read data
    print("Reading TensorBoard data...")
    data = read_tensorboard_logs(args.log_dir)
    
    # Print summary
    print_summary_stats(data)
    
    # Create plots
    plot_path = os.path.join(args.output_dir, "tensorboard_analysis.png")
    fig = plot_training_curves(data, plot_path)
    
    # Export to CSV
    csv_path = os.path.join(args.output_dir, "tensorboard_data.csv")
    df, wide_df = export_to_csv(data, csv_path)
    
    # Show plot if requested
    if args.show_plot:
        plt.show()
    
    print(f"\n✅ Analysis complete! Files saved in: {args.output_dir}")

if __name__ == "__main__":
    main()
