#!/usr/bin/env python3
"""
Extract embeddings from IXI dataset using pretrained MedicalNet ResNet.
This script is specifically designed for the IXI dataset and uses pretrained weights
without training, just for embedding extraction.
"""

import argparse
import torch
import pandas as pd
import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
from pathlib import Path
from train_resnet3d_lightning import BrainMRIModule, BrainMRIDataModule
from monai.data import Dataset
from torch.utils.data import DataLoader
import json
import os


def create_args_for_pretrained(model_name="resnet18", task="regression"):
    """Create args object for pretrained MedicalNet model"""
    class Args:
        pass
    
    args = Args()
    
    # Model configuration
    args.model_name = model_name
    args.pretrained = True
    args.task = task
    args.num_classes = 2  # Not used for regression
    
    # Data configuration
    args.img_col = "path"
    args.y_col = "y"
    args.id_col = "id"
    
    # Head configuration - use simple head for pretrained
    args.simple_head = True
    args.head_dropout = 0.2
    
    # Preprocessing configuration
    args.pixdim = [2.0, 2.0, 2.0]
    args.roi = [128, 128, 128]
    args.no_reorient = False
    args.no_resample = False
    args.no_cropfg = False
    
    # Training configuration (not used but required)
    args.lr = 3e-4
    args.weight_decay = 1e-3
    args.batch_size = 4
    args.epochs = 1
    args.val_frac = 0.1
    args.seed = 1337
    args.num_workers = 0
    args.stratify = False
    args.balanced_sampler = False
    
    # Fine-tuning configuration (not used)
    args.freeze_backbone = False
    args.trainable_layers = "fc"
    args.unfreeze_from = None
    args.lr_head = None
    args.lr_backbone = None
    args.unfreeze_after_epochs = 0
    args.target_norm = "none"  # Don't normalize for embedding extraction
    
    return args


def main():
    parser = argparse.ArgumentParser(description="Extract embeddings from IXI dataset using pretrained MedicalNet")
    parser.add_argument("--input_csv", required=True, help="Path to IXI CSV file (e.g., all_labeled.csv)")
    parser.add_argument("--output_dir", required=True, help="Output directory for embeddings")
    parser.add_argument("--model_name", default="resnet18", 
                       choices=["resnet10", "resnet18", "resnet34", "resnet50", "resnet101", "resnet152", "resnet200"],
                       help="ResNet model architecture")
    parser.add_argument("--batch_size", type=int, default=4, help="Batch size for inference")
    parser.add_argument("--device", default="auto", help="Device to use (cuda/cpu/auto)")
    
    args = parser.parse_args()
    
    # Setup device
    if args.device == "auto":
        device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    else:
        device = torch.device(args.device)
    
    print(f"🧠 Extracting IXI embeddings using pretrained MedicalNet")
    print(f"   Model: {args.model_name}")
    print(f"   Input CSV: {args.input_csv}")
    print(f"   Output dir: {args.output_dir}")
    print(f"   Device: {device}")
    
    # Create model args
    model_args = create_args_for_pretrained(args.model_name, "regression")
    
    # Create model with pretrained weights
    print(f"📥 Loading pretrained MedicalNet {args.model_name}...")
    model = BrainMRIModule(model_args, init_bias=False)
    model = model.to(device)
    model.eval()
    
    print(f"✅ Model loaded successfully")
    
    # Create data module
    data_module = BrainMRIDataModule(model_args)
    
    # Load data
    print(f"📂 Loading data from: {args.input_csv}")
    items, _ = data_module._read_items(
        args.input_csv, model_args.img_col, model_args.y_col, 
        model_args.id_col, model_args.task
    )
    
    dataset = Dataset(items, data_module.transforms)
    dataloader = DataLoader(
        dataset, 
        batch_size=args.batch_size,
        shuffle=False,
        num_workers=0,  # Use 0 workers to avoid issues
        pin_memory=False
    )
    
    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Create output paths
    embeddings_path = output_dir / "embeddings.parquet"
    predictions_path = output_dir / "predictions.parquet"
    index_path = output_dir / "embeddings_index.csv"
    
    print(f"🔄 Extracting embeddings and predictions for {len(items)} samples...")
    print(f"   Using batch_size={args.batch_size}, num_workers=0 for stability")
    
    # Extract embeddings and predictions
    model.export_embeddings_and_predictions(
        dataloader, str(embeddings_path), str(predictions_path), str(index_path)
    )
    
    # Save configuration
    config = {
        'model_name': args.model_name,
        'pretrained': True,
        'task': 'regression',
        'simple_head': True,
        'input_csv': str(args.input_csv),
        'num_samples': len(items),
        'embedding_dim': 512,  # MedicalNet backbone output with simple head
        'device': str(device)
    }
    
    config_path = output_dir / "config.json"
    with open(config_path, "w") as f:
        json.dump(config, f, indent=2)
    
    print(f"✅ Extraction completed!")
    print(f"   Embeddings: {embeddings_path}")
    print(f"   Predictions: {predictions_path}")
    print(f"   Index: {index_path}")
    print(f"   Config: {config_path}")
    
    # Load and print shapes for verification
    embeddings_table = pq.read_table(embeddings_path)
    embeddings = embeddings_table.to_pandas().values
    predictions_table = pq.read_table(predictions_path)
    predictions = predictions_table.to_pandas().values
    index_df = pd.read_csv(index_path)
    
    print(f"\n📊 Results summary:")
    print(f"   Embeddings shape: {embeddings.shape}")
    print(f"   Predictions shape: {predictions.shape}")
    print(f"   Index entries: {len(index_df)}")
    print(f"   Embedding dimension: {embeddings.shape[1] if len(embeddings.shape) > 1 else 'N/A'}")


if __name__ == "__main__":
    main()
