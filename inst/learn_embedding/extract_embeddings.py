#!/usr/bin/env python3
"""
Extract embeddings from brain MRI images using either:
1. A trained model checkpoint, or
2. Pretrained MedicalNet weights

Usage:
  # From trained checkpoint:
  python extract_embeddings.py --checkpoint path/to/best_model.ckpt --input_csv data.csv --output_dir output/
  
  # From trained checkpoint with config name (appends to output_dir):
  python extract_embeddings.py --checkpoint path/to/best_model.ckpt --input_csv data.csv --output_dir output/ --config_name frozen
  
  # From pretrained MedicalNet:
  python extract_embeddings.py --pretrained --model_name resnet18 --input_csv data.csv --output_dir output/
"""

import argparse
import json
import torch
import pandas as pd
import numpy as np
from pathlib import Path
from train_resnet3d_lightning import BrainMRIModule, BrainMRIDataModule
from monai.data import Dataset
from torch.utils.data import DataLoader


def load_from_checkpoint(checkpoint_path, device):
    """Load trained model from checkpoint"""
    checkpoint = torch.load(checkpoint_path, map_location=device, weights_only=False)
    
    # Create args from saved hyperparameters
    class Args:
        pass
    args = Args()
    for key, value in checkpoint['hyper_parameters'].items():
        setattr(args, key, value)
    
    # Create model
    model = BrainMRIModule(args, init_bias=False)
    
    # Register normalization buffers before loading if they exist in checkpoint
    state_dict = checkpoint["state_dict"]
    if "y_mean_buf" in state_dict and "y_std_buf" in state_dict:
        # Register placeholder buffers so load_state_dict can populate them
        model.register_buffer("y_mean_buf", torch.tensor([0.0], dtype=torch.float32))
        model.register_buffer("y_std_buf", torch.tensor([1.0], dtype=torch.float32))
    
    # Load state dict
    missing, unexpected = model.load_state_dict(state_dict, strict=False)
    
    # Print loading info
    print(f"✅ Loaded checkpoint: {checkpoint_path}")
    if hasattr(model, "y_mean_buf") and hasattr(model, "y_std_buf"):
        print(f"   Normalization: y_mean={float(model.y_mean_buf):.3f}, y_std={float(model.y_std_buf):.3f}")
    if missing:
        print(f"   Missing keys: {missing}")
    if unexpected:
        print(f"   Unexpected keys: {unexpected}")
    
    return model, args


def create_pretrained_args(model_name="resnet18", task="regression"):
    """Create args for pretrained MedicalNet model"""
    class Args:
        pass
    args = Args()
    
    # Model config
    args.model_name = model_name
    args.pretrained = True
    args.task = task
    args.num_classes = 2
    args.simple_head = True
    args.head_dropout = 0.2
    
    # Preprocessing config
    args.pixdim = [2.0, 2.0, 2.0]
    args.roi = [128, 128, 128]
    args.no_reorient = False
    args.no_resample = False
    args.no_cropfg = False
    
    # Data config
    args.img_col = "path"
    args.y_col = "y"
    args.id_col = "id"
    
    # Required but unused for extraction
    args.lr = 3e-4
    args.weight_decay = 1e-3
    args.epochs = 1
    args.val_frac = 0.1
    args.seed = 1337
    args.stratify = False
    args.balanced_sampler = False
    args.freeze_backbone = False
    args.trainable_layers = "fc"
    args.unfreeze_from = None
    args.lr_head = None
    args.lr_backbone = None
    args.unfreeze_after_epochs = 0
    args.target_norm = "none"
    
    return args


def main():
    parser = argparse.ArgumentParser(description="Extract embeddings from brain MRI")
    parser.add_argument("--checkpoint", help="Path to trained model checkpoint (.ckpt)")
    parser.add_argument("--pretrained", action="store_true", help="Use pretrained MedicalNet weights")
    parser.add_argument("--model_name", default="resnet18", help="Model architecture (for pretrained)")
    parser.add_argument("--input_csv", required=True, help="CSV with id and path columns")
    parser.add_argument("--output_dir", required=True, help="Output directory")
    parser.add_argument("--config_name", help="Training config descriptor (e.g., 'frozen', 'finetuned', 'scratch') - appended to output_dir")
    parser.add_argument("--batch_size", type=int, default=4, help="Batch size")
    parser.add_argument("--num_workers", type=int, default=0, help="DataLoader workers")
    
    args = parser.parse_args()
    
    # Validate: must specify either checkpoint or pretrained
    if not args.checkpoint and not args.pretrained:
        parser.error("Must specify either --checkpoint or --pretrained")
    if args.checkpoint and args.pretrained:
        parser.error("Cannot specify both --checkpoint and --pretrained")
    
    # Setup output directory with optional config name
    output_dir = Path(args.output_dir)
    if args.config_name:
        output_dir = output_dir / args.config_name
        print(f"📁 Output directory: {output_dir} (with config: {args.config_name})")
    else:
        print(f"📁 Output directory: {output_dir}")
    
    # Create output directory early (needed for temp files)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Setup device
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"🧠 Device: {device}")
    
    # Load model
    if args.checkpoint:
        print(f"📥 Loading from checkpoint...")
        model, model_args = load_from_checkpoint(args.checkpoint, device)
    else:
        print(f"📥 Loading pretrained MedicalNet {args.model_name}...")
        model_args = create_pretrained_args(args.model_name)
        model = BrainMRIModule(model_args, init_bias=False)
        print(f"✅ Pretrained model loaded")
    
    model = model.to(device)
    model.eval()
    
    # Setup data
    data_module = BrainMRIDataModule(model_args)
    
    # Show preprocessing configuration
    print(f"🔧 Preprocessing config:")
    print(f"   Spacing (pixdim): {model_args.pixdim}")
    print(f"   ROI size: {model_args.roi}")
    print(f"   Reorient: {not model_args.no_reorient}")
    print(f"   Resample: {not model_args.no_resample}")
    print(f"   Crop foreground: {not model_args.no_cropfg}")
    
    print(f"📂 Loading data: {args.input_csv}")
    
    # Check if y column exists, if not add dummy column
    df = pd.read_csv(args.input_csv)
    if model_args.y_col not in df.columns:
        print(f"   Note: '{model_args.y_col}' column not found, using dummy values (not needed for extraction)")
        df[model_args.y_col] = 0.0  # dummy value
        temp_csv = output_dir / "temp_input_with_y.csv"
        df.to_csv(temp_csv, index=False)
        input_csv_path = str(temp_csv)
    else:
        input_csv_path = args.input_csv
    
    items, _ = data_module._read_items(
        input_csv_path, model_args.img_col, model_args.y_col,
        model_args.id_col, model_args.task
    )
    
    dataset = Dataset(items, data_module.transforms)
    dataloader = DataLoader(
        dataset, batch_size=args.batch_size, shuffle=False,
        num_workers=args.num_workers, pin_memory=False
    )
    
    # Setup output paths
    embeddings_path = output_dir / "embeddings.npy"
    predictions_path = output_dir / "predictions.npy"
    index_path = output_dir / "embeddings_index.csv"
    config_path = output_dir / "extraction_config.json"
    
    # Extract embeddings
    print(f"🔄 Extracting embeddings for {len(items)} samples...")
    model.export_embeddings_and_predictions(
        dataloader, str(embeddings_path), str(predictions_path), str(index_path)
    )
    
    # Save extraction configuration
    extraction_config = {
        'input_csv': str(args.input_csv),
        'num_samples': len(items),
        'batch_size': args.batch_size,
        'config_name': args.config_name,
        'model_source': 'checkpoint' if args.checkpoint else 'pretrained',
    }
    
    if args.checkpoint:
        extraction_config['checkpoint_path'] = str(args.checkpoint)
        # Add training configuration info from loaded model
        if hasattr(model_args, 'freeze_backbone'):
            extraction_config['freeze_backbone'] = model_args.freeze_backbone
        if hasattr(model_args, 'trainable_layers'):
            extraction_config['trainable_layers'] = model_args.trainable_layers
        if hasattr(model_args, 'unfreeze_from'):
            extraction_config['unfreeze_from'] = model_args.unfreeze_from
        if hasattr(model_args, 'pretrained'):
            extraction_config['pretrained_backbone'] = model_args.pretrained
        if hasattr(model_args, 'simple_head'):
            extraction_config['simple_head'] = model_args.simple_head
    else:
        extraction_config['model_name'] = args.model_name
        extraction_config['pretrained'] = True
    
    with open(config_path, 'w') as f:
        json.dump(extraction_config, f, indent=2)
    
    # Clean up temporary CSV if created
    temp_csv = output_dir / "temp_input_with_y.csv"
    if temp_csv.exists():
        temp_csv.unlink()
    
    # Print summary
    embeddings = np.load(embeddings_path)
    predictions = np.load(predictions_path)
    
    print(f"\n✅ Done!")
    print(f"   Embeddings: {embeddings_path} (shape: {embeddings.shape})")
    print(f"   Predictions: {predictions_path} (shape: {predictions.shape})")
    print(f"   Index: {index_path}")
    print(f"   Config: {config_path}")


if __name__ == "__main__":
    main()

