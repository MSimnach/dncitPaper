#!/usr/bin/env python3
"""
Simple wrapper script to run IXI embedding extraction that can be called from R.
This avoids complex import issues by keeping everything in one file.
"""

import sys
import os

def run_ixi_extraction(input_csv, output_dir, model_name="resnet18", batch_size=4):
    """
    Run IXI embedding extraction with the given parameters.
    """
    # Add the current directory to Python path
    current_dir = os.path.dirname(os.path.abspath(__file__))
    if current_dir not in sys.path:
        sys.path.insert(0, current_dir)
    
    try:
        # Import required modules
        import torch
        
        # Optimize CPU usage
        torch.set_num_threads(min(16, os.cpu_count() or 1))  # Limit CPU threads to avoid oversubscription
        
        print(f"Using PyTorch {torch.__version__}")
        print(f"CUDA available: {torch.cuda.is_available()}")
        print(f"CPU threads: {torch.get_num_threads()}")
        
        # Import the extraction script
        from extract_ixi_embeddings import create_args_for_pretrained
        from train_resnet3d_lightning import BrainMRIModule, BrainMRIDataModule
        
        # Setup device
        device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        
        print(f'🧠 Extracting IXI embeddings using pretrained MedicalNet')
        print(f'   Model: {model_name}')
        print(f'   Input CSV: {input_csv}')
        print(f'   Output dir: {output_dir}')
        print(f'   Device: {device}')
        
        # Create model args
        model_args = create_args_for_pretrained(model_name, 'regression')
        
        # Create model with pretrained weights
        print(f'📥 Loading pretrained MedicalNet {model_name}...')
        model = BrainMRIModule(model_args, init_bias=False)
        model = model.to(device)
        model.eval()
        
        # Optimize model for inference
        use_half_precision = False  # Disable FP16 for now due to compatibility issues
        if torch.cuda.is_available() and use_half_precision:
            # Enable mixed precision for faster inference on GPU
            model = model.half()  # Use FP16 for faster inference
            print(f'⚡ Enabled FP16 (half precision) for faster GPU inference')
        
        # Disable gradient computation for inference
        torch.set_grad_enabled(False)
        
        
        print(f'✅ Model loaded and optimized for inference')
        
        # Create data module
        data_module = BrainMRIDataModule(model_args)
        
        # Load data
        print(f'📂 Loading data from: {input_csv}')
        items, _ = data_module._read_items(
            input_csv, model_args.img_col, model_args.y_col, 
            model_args.id_col, model_args.task
        )
        
        from monai.data import Dataset
        from torch.utils.data import DataLoader
        
        dataset = Dataset(items, data_module.transforms)
        # Optimize DataLoader for faster inference
        num_workers = min(8, os.cpu_count() or 1)  # Use multiple workers but not too many
        pin_memory = torch.cuda.is_available()  # Enable pin_memory if CUDA available
        
        dataloader = DataLoader(
            dataset, 
            batch_size=batch_size,
            shuffle=False,
            num_workers=num_workers,
            pin_memory=pin_memory,
            persistent_workers=num_workers > 0,  # Keep workers alive between epochs
            prefetch_factor=2 if num_workers > 0 else 2  # Prefetch batches
        )
        
        print(f'⚡ DataLoader optimizations: batch_size={batch_size}, num_workers={num_workers}, pin_memory={pin_memory}')
        
        # Create output directory
        from pathlib import Path
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Create output paths
        embeddings_path = output_dir / 'embeddings.parquet'
        predictions_path = output_dir / 'predictions.parquet'
        index_path = output_dir / 'embeddings_index.csv'
        
        print(f'🔄 Extracting embeddings and predictions for {len(items)} samples...')
        print(f'⚡ Using optimized settings: batch_size={batch_size}, workers={num_workers}')
        
        # Extract embeddings and predictions with error handling
        import time
        start_time = time.time()
        
        try:
            model.export_embeddings_and_predictions(
                dataloader, str(embeddings_path), str(predictions_path), str(index_path)
            )
        except RuntimeError as e:
            if "out of memory" in str(e).lower() or "cuda" in str(e).lower():
                print(f"⚠️  Memory error with batch_size={batch_size}, trying smaller batch...")
                # Recreate dataloader with smaller batch size
                smaller_batch = max(1, batch_size // 2)
                dataloader = DataLoader(
                    dataset, 
                    batch_size=smaller_batch,
                    shuffle=False,
                    num_workers=min(4, num_workers),  # Also reduce workers
                    pin_memory=pin_memory,
                    persistent_workers=min(4, num_workers) > 0,
                    prefetch_factor=2 if min(4, num_workers) > 0 else 2
                )
                print(f"🔄 Retrying with batch_size={smaller_batch}")
                model.export_embeddings_and_predictions(
                    dataloader, str(embeddings_path), str(predictions_path), str(index_path)
                )
            else:
                raise e
        
        elapsed_time = time.time() - start_time
        samples_per_sec = len(items) / elapsed_time if elapsed_time > 0 else 0
        print(f'⏱️  Processing completed in {elapsed_time:.1f}s ({samples_per_sec:.1f} samples/sec)')
        
        # Save configuration
        import json
        config = {
            'model_name': model_name,
            'pretrained': True,
            'task': 'regression',
            'simple_head': True,
            'input_csv': str(input_csv),
            'num_samples': len(items),
            'embedding_dim': 512,
            'device': str(device)
        }
        
        config_path = output_dir / 'config.json'
        with open(config_path, 'w') as f:
            json.dump(config, f, indent=2)
        
        print(f'✅ Extraction completed!')
        print(f'   Embeddings: {embeddings_path}')
        print(f'   Predictions: {predictions_path}')
        print(f'   Index: {index_path}')
        
        return True
        
    except Exception as e:
        print(f"❌ Error during extraction: {str(e)}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == "__main__":
    # For command line usage
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_csv", required=True)
    parser.add_argument("--output_dir", required=True)
    parser.add_argument("--model_name", default="resnet18")
    parser.add_argument("--batch_size", type=int, default=4)
    
    args = parser.parse_args()
    
    success = run_ixi_extraction(
        args.input_csv, 
        args.output_dir, 
        args.model_name, 
        args.batch_size
    )
    
    sys.exit(0 if success else 1)
