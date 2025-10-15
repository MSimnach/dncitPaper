#!/usr/bin/env python3
"""
Simple script to extract TensorBoard data to pandas DataFrame
"""

import os
import pandas as pd
import numpy as np
from tensorboard.backend.event_processing.event_accumulator import EventAccumulator

def read_tensorboard_simple(log_dir):
    """Simple function to read TensorBoard data into pandas DataFrame"""
    
    # Find latest version
    versions = [d for d in os.listdir(log_dir) if d.startswith('version_')]
    latest_version = sorted(versions, key=lambda x: int(x.split('_')[1]))[-1]
    event_dir = os.path.join(log_dir, latest_version)
    
    # Load events
    event_acc = EventAccumulator(event_dir)
    event_acc.Reload()
    
    # Extract epoch-level metrics (most important)
    epoch_data = {}
    epoch_metrics = ['train_loss_epoch', 'val_loss_epoch', 'loss_diff']
    
    for metric in epoch_metrics:
        if metric in event_acc.Tags()['scalars']:
            scalar_events = event_acc.Scalars(metric)
            steps = [event.step for event in scalar_events]
            values = [event.value for event in scalar_events]
            epoch_data[metric] = values
            if 'epoch' not in epoch_data:
                epoch_data['epoch'] = steps
    
    return pd.DataFrame(epoch_data)

# Example usage
if __name__ == "__main__":
    import sys
    if len(sys.argv) < 2:
        print("Usage: python simple_tensorboard_reader.py <tensorboard_log_dir>")
        sys.exit(1)
    
    log_dir = sys.argv[1]
    df = read_tensorboard_simple(log_dir)
    
    print("Epoch-level training data:")
    print(df)
    print(f"\nData shape: {df.shape}")
    print(f"Columns: {list(df.columns)}")
    
    # Save to CSV
    output_file = "training_epochs.csv"
    df.to_csv(output_file, index=False)
    print(f"\nData saved to: {output_file}")
