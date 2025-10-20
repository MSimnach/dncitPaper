#!/usr/bin/env python3
"""
Complete pipeline for training ResNet on brain MRIs with train/test split.

This script:
1. Splits the data into train/test (if not already done)
2. Trains ResNet on train set (with train/val split) with optional MedicalNet pretrained weights
   - Supports fine-tuning with layer freezing, discriminative learning rates, and gradual unfreezing
3. Extracts penultimate layer embeddings and predictions from test set
"""

import argparse
import subprocess
import sys
from pathlib import Path


def run_command(cmd, description, env=None):
    """Run a shell command and handle errors"""
    print(f"\n🔄 {description}")
    print(f"Command: {' '.join(cmd)}")
    
    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True, env=env)
        print(f"✅ {description} completed successfully")
        if result.stdout:
            print("Output:", result.stdout)
        return True
    except subprocess.CalledProcessError as e:
        print(f"❌ {description} failed with exit code {e.returncode}")
        if e.stdout:
            print("STDOUT:", e.stdout)
        if e.stderr:
            print("STDERR:", e.stderr)
        return False


def main():
    parser = argparse.ArgumentParser(description="Complete train/test pipeline for ResNet brain MRI")
    
    # Data paths
    parser.add_argument("--input_csv", required=True,
                       help="Path to all_labeled.csv file")
    parser.add_argument("--output_dir", required=True,
                       help="Output directory for training results")
    parser.add_argument("--script_dir", required=True,
                       help="Directory containing the split_data.py and train_resnet3d_lightning.py scripts")
    parser.add_argument("--id_col", default=None,
                       help="Column name for subject IDs (if not provided, uses file path stem)")
    
    # Data splitting
    parser.add_argument("--test_size", type=float, default=0.5,
                       help="Fraction of data for test set (default: 0.5)")
    parser.add_argument("--skip_split", action="store_true",
                       help="Skip data splitting (assume train/test CSVs already exist)")
    
    # Training parameters
    parser.add_argument("--model_name", default="resnet18",
                       choices=["resnet10", "resnet18", "resnet34", "resnet50", "resnet101", "resnet152", "resnet200"])
    parser.add_argument("--pretrained", action="store_true",
                       help="Use MedicalNet pretrained weights")
    parser.add_argument("--task", choices=["regression", "binary", "multiclass"], default="regression")
    parser.add_argument("--epochs", type=int, default=30)
    parser.add_argument("--batch_size", type=int, default=6)
    parser.add_argument("--lr", type=float, default=3e-4)
    parser.add_argument("--val_frac", type=float, default=0.2,
                       help="Validation fraction from training set")
    
    # Fine-tuning parameters
    parser.add_argument("--freeze_backbone", action="store_true",
                       help="Freeze all backbone layers (conv1..layer4); train head only")
    parser.add_argument("--simple_head", action="store_true",
                       help="Use simple Linear head instead of MLP")
    parser.add_argument("--trainable_layers", type=str, default="fc",
                       help="Comma-separated list of layers to unfreeze: e.g. 'fc' or 'layer4,fc' or 'layer3,layer4,fc'.")
    parser.add_argument("--unfreeze_from", type=str, default=None,
                       help="Alternative to trainable_layers: unfreeze this block and all later blocks (one of: conv1,layer1,layer2,layer3,layer4,fc).")
    parser.add_argument("--lr_head", type=float, default=None,
                       help="Learning rate for head (fc). Defaults to --lr if not set.")
    parser.add_argument("--lr_backbone", type=float, default=None,
                       help="Learning rate for unfrozen backbone blocks. Defaults to 0.1 * --lr if not set.")
    parser.add_argument("--head_dropout", type=float, default=0.2,
                       help="Dropout rate for MLP head")
    parser.add_argument("--unfreeze_after_epochs", type=int, default=0,
                       help="Gradual unfreezing: after this many epochs, also unfreeze the next earlier block.")
    
    # System
    parser.add_argument("--num_workers", type=int, default=8)
    parser.add_argument("--seed", type=int, default=1337)
    parser.add_argument("--amp", action="store_true", help="Use mixed precision training")
    
    # Logging
    parser.add_argument("--use_wandb", action="store_true")
    parser.add_argument("--use_tensorboard", action="store_true")
    
    args = parser.parse_args()
    
    # Convert paths
    input_csv = Path(args.input_csv)
    output_dir = Path(args.output_dir)
    csv_dir = input_csv.parent
    script_dir = Path(args.script_dir)
    
    # Create output directory
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print("🧠 ResNet Brain MRI Train/Test Pipeline")
    print(f"   Input CSV: {input_csv}")
    print(f"   Output dir: {output_dir}")
    print(f"   Model: {args.model_name}")
    print(f"   Pretrained: {'MedicalNet' if args.pretrained else 'Random init'}")
    print(f"   Task: {args.task}")
    
    # Step 1: Split data (if needed)
    train_csv = csv_dir / "train_labeled.csv"
    test_csv = csv_dir / "test_labeled.csv"
    
    if not args.skip_split or not (train_csv.exists() and test_csv.exists()):
        split_cmd = [
            sys.executable, str(script_dir / "split_data.py"),
            "--input_csv", str(input_csv),
            "--output_dir", str(csv_dir),
            "--test_size", str(args.test_size),
            "--seed", str(args.seed)
        ]
        
        if not run_command(split_cmd, "Data splitting"):
            return 1
    else:
        print(f"✅ Using existing train/test splits:")
        print(f"   Train: {train_csv}")
        print(f"   Test: {test_csv}")
    
    # Step 2: Train model with embedding extraction
    # Set environment to use GPU 1 (GPU 0 is often occupied)
    import os
    #os.environ["CUDA_VISIBLE_DEVICES"] = "0"
    
    train_cmd = [
        sys.executable, str(script_dir / "train_resnet3d_lightning.py"),
        "--uaux_csv", str(input_csv),  # This will be overridden by --use_train_test_split
        "--use_train_test_split",
        "--out_dir", str(output_dir),
        "--model_name", args.model_name,
        "--task", args.task,
        "--epochs", str(args.epochs),
        "--batch_size", str(args.batch_size),
        "--lr", str(args.lr),
        "--val_frac", str(args.val_frac),
        "--num_workers", str(args.num_workers),
        "--seed", str(args.seed)
    ]
    
    # Add optional arguments
    if args.pretrained:
        train_cmd.append("--pretrained")
    if args.id_col:
        train_cmd.extend(["--id_col", args.id_col])
    if args.freeze_backbone:
        train_cmd.append("--freeze_backbone")
    if args.simple_head:
        train_cmd.append("--simple_head")
    if args.lr_head is not None:
        train_cmd.extend(["--lr_head", str(args.lr_head)])
    if args.lr_backbone is not None:
        train_cmd.extend(["--lr_backbone", str(args.lr_backbone)])
    if args.head_dropout != 0.2:  # Only add if different from default
        train_cmd.extend(["--head_dropout", str(args.head_dropout)])
    if args.trainable_layers != "fc":  # Only add if different from default
        train_cmd.extend(["--trainable_layers", args.trainable_layers])
    if args.unfreeze_from:
        train_cmd.extend(["--unfreeze_from", args.unfreeze_from])
    if args.unfreeze_after_epochs > 0:
        train_cmd.extend(["--unfreeze_after_epochs", str(args.unfreeze_after_epochs)])
    
    # Set environment variable for the subprocess
    env = os.environ.copy()
    #env["CUDA_VISIBLE_DEVICES"] = "0"
    
    # Add optional flags
    if args.amp:
        train_cmd.append("--amp")
    if args.use_wandb:
        train_cmd.append("--use_wandb")
    if args.use_tensorboard:
        train_cmd.append("--use_tensorboard")
    
    if not run_command(train_cmd, "Model training and embedding extraction", env=env):
        return 1
    
    print(f"\n🎉 Pipeline completed successfully!")
    print(f"   Training results: {output_dir}")
    print(f"   Test embeddings: {output_dir / 'test_embeddings.npy'}")
    print(f"   Test predictions: {output_dir / 'test_predictions.npy'}")
    print(f"   Embedding index: {output_dir / 'test_embeddings_index.csv'}")
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
