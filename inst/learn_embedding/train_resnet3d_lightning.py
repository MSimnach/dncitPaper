#!/usr/bin/env python3
# Train a 3D ResNet (MONAI) from scratch on U_aux to predict Y from X (T1 MRI).
# - Supports regression or (binary/multiclass) classification
# - PyTorch Lightning for clean training and automatic progress visualization
# - Early stopping, AMP, and optional export of penultimate embeddings

import os, csv, json, argparse, math, random, time
from pathlib import Path
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import WeightedRandomSampler, DataLoader
import pytorch_lightning as pl
from pytorch_lightning.callbacks import EarlyStopping, ModelCheckpoint, RichProgressBar
from pytorch_lightning.loggers import CSVLogger

from monai.utils import set_determinism
from monai.data import Dataset
from monai.transforms import (
    Compose, LoadImaged, EnsureChannelFirstd, Orientationd, Spacingd,
    CropForegroundd, ResizeWithPadOrCropd, ScaleIntensityRangePercentilesd, ToTensord
)
from monai.networks.nets import (
    resnet10, resnet18, resnet34, resnet50, resnet101, resnet152, resnet200
)
from sklearn.model_selection import train_test_split

# Optional loggers
try:
    from pytorch_lightning.loggers import WandbLogger
    HAS_WANDB = True
except ImportError:
    HAS_WANDB = False

try:
    from pytorch_lightning.loggers import TensorBoardLogger
    HAS_TENSORBOARD = True
except ImportError:
    HAS_TENSORBOARD = False

MODEL_FNS = {
    "resnet10": resnet10, "resnet18": resnet18, "resnet34": resnet34,
    "resnet50": resnet50, "resnet101": resnet101, "resnet152": resnet152, "resnet200": resnet200,
}

# Order of backbone layers for freeze/unfreeze operations
BACKBONE_ORDER = ["conv1", "bn1", "act", "maxpool", "layer1", "layer2", "layer3", "layer4", "avgpool", "fc"]

class BrainMRIModule(pl.LightningModule):
    """PyTorch Lightning module for 3D ResNet brain MRI training"""
    
    def __init__(self, args, init_bias=True):
        super().__init__()
        self.args = args
        self.save_hyperparameters(vars(args))
        
        # Build model
        self.model = self._build_model()
        self._apply_freeze_plan()
        self.criterion = self._get_criterion()
        
        # Log trainable layers info
        self._log_trainable_layers()
        
        # Initialize bias only if requested (not when loading from checkpoint)
        self._init_bias = init_bias
        
        # Track metrics for logging
        self.train_losses = []
        self.val_losses = []
        
    def _build_model(self):
        if self.args.task == "regression":
            out_ch = 1
        elif self.args.task == "binary":
            out_ch = 1
        else:
            out_ch = self.args.num_classes

        fn = MODEL_FNS[self.args.model_name]
        
        # MedicalNet pretrained weights require specific parameters
        if self.args.pretrained:
            model = fn(
                spatial_dims=3,
                n_input_channels=1,
                num_classes=out_ch,
                pretrained=self.args.pretrained,
                feed_forward=False,  # Required for MedicalNet pretrained weights
                shortcut_type="A",   # Required for MedicalNet pretrained weights
                bias_downsample=True,  # Required for MedicalNet pretrained weights
            )
            
            # MedicalNet provides a 512-dimensional embedding from its backbone
            feature_dim = 512  # MedicalNet backbone output dimension
            
            if self.args.simple_head:
                # Simple Linear head: 512 -> output
                model.fc = torch.nn.Linear(feature_dim, out_ch)
                print(f"[debug] Added simple head: Linear({feature_dim} -> {out_ch})")
            else:
                # 2-layer MLP with configurable dropout: 512 -> 512 -> 256 -> output
                model.fc = torch.nn.Sequential(
                    torch.nn.Linear(feature_dim, 512),
                    torch.nn.ReLU(),
                    torch.nn.Dropout(self.args.head_dropout),
                    torch.nn.Linear(512, 256),
                    torch.nn.ReLU(),
                    torch.nn.Dropout(self.args.head_dropout),
                    torch.nn.Linear(256, out_ch)
                )
                print(f"[debug] Added 2-layer MLP: {feature_dim} -> 512 -> 256 -> {out_ch} (dropout={self.args.head_dropout})")
            
            return model
        else:
            return fn(
                spatial_dims=3,
                n_input_channels=1,
                num_classes=out_ch,
                pretrained=False,
            )
    
    def _set_requires_grad(self, module, requires_grad: bool):
        """Set requires_grad for all parameters in module and handle BatchNorm eval mode"""
        for p in module.parameters():
            p.requires_grad = requires_grad
        # If freezing, put BatchNorms in eval so running stats are not updated
        if not requires_grad:
            for m in module.modules():
                if isinstance(m, (nn.BatchNorm3d, nn.InstanceNorm3d)):
                    m.eval()

    def _freeze_all(self):
        """Freeze all backbone layers"""
        for name in BACKBONE_ORDER:
            if hasattr(self.model, name):
                self._set_requires_grad(getattr(self.model, name), False)

    def _unfreeze_layers(self, layers):
        """Unfreeze specific layers"""
        for name in layers:
            if hasattr(self.model, name):
                self._set_requires_grad(getattr(self.model, name), True)

    def _apply_freeze_plan(self):
        """Apply freeze/unfreeze based on args."""
        if not self.args.freeze_backbone and self.args.unfreeze_from is None and self.args.trainable_layers == "fc":
            # default behavior: train everything if not freezing; else fall through
            return

        # Start: freeze everything
        self._freeze_all()

        # Strategy A: unfreeze specific layers list
        if self.args.trainable_layers and self.args.unfreeze_from is None:
            layers = [s.strip() for s in self.args.trainable_layers.split(",") if s.strip()]
            self._unfreeze_layers(layers)
            return

        # Strategy B: unfreeze from a block onward
        if self.args.unfreeze_from is not None:
            if self.args.unfreeze_from not in BACKBONE_ORDER:
                raise ValueError(f"--unfreeze_from must be one of {BACKBONE_ORDER}")
            start_idx = BACKBONE_ORDER.index(self.args.unfreeze_from)
            layers = BACKBONE_ORDER[start_idx:]
            self._unfreeze_layers(layers)
    
    def _log_trainable_layers(self):
        """Log information about which layers are trainable"""
        trainable_layers = []
        frozen_layers = []
        
        for name in BACKBONE_ORDER:
            if hasattr(self.model, name):
                module = getattr(self.model, name)
                if any(p.requires_grad for p in module.parameters()):
                    trainable_layers.append(name)
                else:
                    frozen_layers.append(name)
        
        if self.args.freeze_backbone or self.args.unfreeze_from or self.args.trainable_layers != "fc":
            print(f"🔧 Fine-tuning setup:")
            if trainable_layers:
                print(f"   ✅ Trainable layers: {', '.join(trainable_layers)}")
            if frozen_layers:
                print(f"   ❄️  Frozen layers: {', '.join(frozen_layers)}")
    
    def _init_head_bias(self):
        """Initialize final layer bias to training mean"""
        try:
            if hasattr(self, 'y_mean_buf') and self.y_mean_buf is not None:
                y_mean = float(self.y_mean_buf.item())
            else:
                y_mean = 50.0  # fallback
        except Exception:
            y_mean = 50.0  # fallback
        
        with torch.no_grad():
            if hasattr(self.model.fc, 'bias') and self.model.fc.bias is not None:
                # Simple linear head
                self.model.fc.bias.fill_(y_mean)
                print(f"🎯 Initialized head bias to {y_mean:.1f}")
            elif hasattr(self.model.fc, '__getitem__'):
                # Sequential MLP - get last layer
                last_layer = self.model.fc[-1]
                if hasattr(last_layer, 'bias') and last_layer.bias is not None:
                    last_layer.bias.fill_(y_mean)
                    # Optionally scale down the last weight matrix
                    last_layer.weight.mul_(0.01)
                    print(f"🎯 Initialized final layer bias to {y_mean:.1f} and scaled weights by 0.01")
    
    def _get_criterion(self):
        if self.args.task == "regression":
            return torch.nn.MSELoss()
        elif self.args.task == "binary":
            return torch.nn.BCEWithLogitsLoss()
        else:
            return torch.nn.CrossEntropyLoss()
    
    def forward(self, x):
        return self.model(x)
    
    def _step(self, batch, batch_idx, stage):
        x = batch["image"]
        y = batch["y"]
        
        # Prepare targets
        if self.args.task == "regression":
            y = torch.as_tensor(y, dtype=torch.float32, device=self.device).view(-1, 1)
            # Apply z-score normalization if enabled
            if (self.args.target_norm == "zscore" and 
                hasattr(self, 'y_mean_buf') and hasattr(self, 'y_std_buf') and
                self.y_mean_buf is not None and self.y_std_buf is not None):
                y = (y - self.y_mean_buf) / (self.y_std_buf + 1e-6)
        elif self.args.task == "binary":
            y = torch.as_tensor(y, dtype=torch.float32, device=self.device).view(-1, 1)
        else:
            y = torch.as_tensor(y, dtype=torch.long, device=self.device)
        
        # Forward pass
        logits = self(x)
        loss = self.criterion(logits, y)
        
        # Log metrics
        self.log(f'{stage}_loss', loss, on_step=True, on_epoch=True, prog_bar=True, logger=True)
        
        return {'loss': loss, 'logits': logits, 'targets': y}
    
    def training_step(self, batch, batch_idx):
        return self._step(batch, batch_idx, 'train')
    
    def validation_step(self, batch, batch_idx):
        return self._step(batch, batch_idx, 'val')
    
    def on_train_epoch_end(self):
        # Get average loss for the epoch
        train_loss = self.trainer.callback_metrics.get('train_loss_epoch', 0)
        self.train_losses.append(float(train_loss))
        
        # Gradual unfreeze: after k epochs, unfreeze the next earlier block
        k = self.args.unfreeze_after_epochs
        if k and self.current_epoch + 1 in [k, 2*k, 3*k, 4*k]:  # stepwise at k,2k,3k,4k
            # Determine which block to unfreeze next (walk backward)
            # Start from the earliest currently-frozen among ['layer3','layer2','layer1','conv1']
            for candidate in ["layer3", "layer2", "layer1", "conv1"]:
                mod = getattr(self.model, candidate, None)
                if mod is None:
                    continue
                # Check if currently frozen (any param)
                frozen = all(not p.requires_grad for p in mod.parameters())
                if frozen:
                    print(f"🔓 Gradual unfreeze: enabling {candidate} at epoch {self.current_epoch+1}")
                    self._set_requires_grad(mod, True)
                    # Rebuild optimizer with new param groups (Lightning will pick it up next step)
                    self.trainer.strategy.optimizer_zero_grad(self.trainer.optimizers[0], self.current_epoch, 0)
                    self.trainer.optimizers = [self.configure_optimizers()]  # naive refresh
                    break
    
    def on_validation_epoch_end(self):
        # Get average loss for the epoch
        val_loss = self.trainer.callback_metrics.get('val_loss_epoch', 0)
        self.val_losses.append(float(val_loss))
        
        # Log additional metrics
        if len(self.train_losses) > 0 and len(self.val_losses) > 0:
            loss_diff = float(val_loss) - self.train_losses[-1]
            self.log('loss_diff', loss_diff, prog_bar=True, logger=True)
    
    def configure_optimizers(self):
        # Check if we're using fine-tuning features (discriminative LRs)
        using_fine_tuning = (
            self.args.freeze_backbone or 
            self.args.unfreeze_from is not None or 
            self.args.trainable_layers != "fc" or
            self.args.lr_head is not None or
            self.args.lr_backbone is not None
        )
        
        if not using_fine_tuning:
            # Default behavior: single learning rate for all parameters (backward compatibility)
            optimizer = torch.optim.AdamW(
                self.parameters(), 
                lr=self.args.lr, 
                weight_decay=self.args.weight_decay
            )
            return optimizer
        
        # Fine-tuning mode: collect params into head and backbone groups
        head = []
        backbone = []
        for name, p in self.model.named_parameters():
            if not p.requires_grad:
                continue
            if name.startswith("fc.") or name == "fc.weight" or name == "fc.bias":
                head.append(p)
            else:
                backbone.append(p)

        lr_head = self.args.lr if self.args.lr_head is None else self.args.lr_head
        lr_back = (0.1 * self.args.lr) if self.args.lr_backbone is None else self.args.lr_backbone

        param_groups = []
        if backbone:
            param_groups.append({"params": backbone, "lr": lr_back})
        if head:
            param_groups.append({"params": head, "lr": lr_head})

        optimizer = torch.optim.AdamW(param_groups, weight_decay=self.args.weight_decay)
        
        # Debug: Print optimizer param groups
        for i, g in enumerate(optimizer.param_groups):
            print(f"[opt] group{i} lr={g['lr']:.2e} params={sum(p.numel() for p in g['params']):,}")
        
        return optimizer
    
    def export_embeddings_and_predictions(self, dataloader, embeddings_path, predictions_path, index_path):
        """Export penultimate layer embeddings and predictions"""
        self.eval()
        feats, preds, sids = [], [], []
        
        with torch.no_grad():
            for batch in dataloader:
                x = batch["image"].to(self.device)
                
                if self.args.pretrained:
                    # For pretrained MedicalNet models, extract features from backbone
                    z = self.model.conv1(x)
                    z = self.model.bn1(z)
                    z = self.model.act(z)  # MONAI uses 'act' instead of 'relu'
                    z = self.model.maxpool(z)
                    z = self.model.layer1(z)
                    z = self.model.layer2(z)
                    z = self.model.layer3(z)
                    z = self.model.layer4(z)
                    gap = self.model.avgpool(z).flatten(1)  # (B, 512)
                    
                    if not self.args.simple_head:
                        # For MLP head, extract features from penultimate MLP layer (256-dim)
                        # self.model.fc is Sequential: Linear(512, 512) -> ReLU -> Dropout -> Linear(512, 256) -> ReLU -> Dropout -> Linear(256, out_ch)
                        mlp_features = self.model.fc[0](gap)  # Linear(512, 512)
                        mlp_features = self.model.fc[1](mlp_features)  # ReLU
                        mlp_features = self.model.fc[2](mlp_features)  # Dropout
                        mlp_features = self.model.fc[3](mlp_features)  # Linear(512, 256)
                        mlp_features = self.model.fc[4](mlp_features)  # ReLU
                        # Stop before final dropout and linear layer to get 256-dim embeddings
                        gap = mlp_features  # (B, 256)
                    # For simple head, keep the 512-dim backbone features as embeddings
                else:
                    # For non-pretrained models, extract features from backbone as before
                    z = self.model.conv1(x)
                    z = self.model.bn1(z)
                    z = self.model.act(z)  # MONAI uses 'act' instead of 'relu'
                    z = self.model.maxpool(z)
                    z = self.model.layer1(z)
                    z = self.model.layer2(z)
                    z = self.model.layer3(z)
                    z = self.model.layer4(z)
                    # Use the model's built-in avgpool instead of F.adaptive_avg_pool3d
                    gap = self.model.avgpool(z).flatten(1)  # (B, C)
                
                # Get predictions - always use forward pass through the PyTorch Lightning module
                # This ensures we get the final output regardless of the underlying model architecture
                logits = self(x)
                
                # Convert logits to predictions based on task
                if self.args.task == "regression":
                    pred = logits  # Keep as is for regression
                    # Denormalize if z-score was used during training
                    if (self.args.target_norm == "zscore" and 
                        hasattr(self, 'y_mean_buf') and hasattr(self, 'y_std_buf') and
                        self.y_mean_buf is not None and self.y_std_buf is not None):
                        pred = pred * (self.y_std_buf + 1e-6) + self.y_mean_buf
                elif self.args.task == "binary":
                    pred = torch.sigmoid(logits)  # Convert to probabilities
                else:  # multiclass
                    pred = F.softmax(logits, dim=1)  # Convert to probabilities
                
                feats.append(gap.cpu().numpy())
                preds.append(pred.cpu().numpy())
                sids.extend(batch["subject_id"])
        
        # Save embeddings, predictions and index
        feats = np.concatenate(feats, 0)
        preds = np.concatenate(preds, 0)
        
        
        Path(embeddings_path).parent.mkdir(parents=True, exist_ok=True)
        np.save(embeddings_path, feats)
        np.save(predictions_path, preds)
        
        with open(index_path, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["row", "subject_id"])
            for i, sid in enumerate(sids):
                writer.writerow([i, sid])
        
        print(f"[export] Saved embeddings: {embeddings_path} shape={feats.shape}")
        print(f"[export] Saved predictions: {predictions_path} shape={preds.shape}")
        print(f"[export] Saved index: {index_path}")
    
    def export_embeddings(self, dataloader, output_path, index_path):
        """Export penultimate layer embeddings (backward compatibility)"""
        # For backward compatibility, extract predictions path from embeddings path
        predictions_path = str(Path(output_path).with_name(Path(output_path).stem + "_predictions.npy"))
        self.export_embeddings_and_predictions(dataloader, output_path, predictions_path, index_path)

class BrainMRIDataModule(pl.LightningDataModule):
    """PyTorch Lightning data module for brain MRI data"""
    
    def __init__(self, args):
        super().__init__()
        self.args = args
        self.transforms = self._build_transforms()
        
    def _build_transforms(self):
        t = [LoadImaged(keys=["image"]), EnsureChannelFirstd(keys=["image"])]
        if not self.args.no_reorient:
            t += [Orientationd(keys=["image"], axcodes="RAS")]
        if not self.args.no_resample:
            t += [Spacingd(keys=["image"], pixdim=tuple(self.args.pixdim), mode=("bilinear"))]
        if not self.args.no_cropfg:
            t += [CropForegroundd(keys=["image"], source_key="image")]
        t += [ResizeWithPadOrCropd(keys=["image"], spatial_size=tuple(self.args.roi))]
        t += [ScaleIntensityRangePercentilesd(
                keys=["image"], lower=1, upper=99, b_min=0.0, b_max=1.0, clip=True
             )]
        t += [ToTensord(keys=["image"])]
        return Compose(t)
    
    def _read_items(self, csv_path, img_col, y_col, id_col=None, task="regression"):
        df = pd.read_csv(csv_path)
        if img_col not in df.columns or y_col not in df.columns:
            raise ValueError(f"CSV must have columns '{img_col}' and '{y_col}'. Found: {df.columns.tolist()}")
        
        items = []
        ys = []
        for _, r in df.iterrows():
            p = str(r[img_col])
            sid = str(r[id_col]) if id_col and id_col in df.columns else Path(p).stem
            y = r[y_col]
            if len(items) < 5:  # Only print for first 5 samples to avoid spam
                print(f"[DEBUG] Sample {len(items)}:")
                print(f"  Image path: {p}")
                print(f"  Subject ID: {sid}")
                print(f"  Raw y value: {y} (type: {type(y)})")
                print(f"  ID column exists: {id_col in df.columns if id_col else 'No id_col specified'}")
            if task == "regression":
                y = float(y)
            elif task == "binary":
                y = int(y)  # 0/1
            else:
                y = int(y)  # 0..K-1
            
            items.append({"image": p, "y": y, "subject_id": sid})
            ys.append(y)
        
        # ADD SUMMARY DEBUG CODE HERE (before return)
        print(f"[DEBUG] Data loading summary:")
        print(f"  Total samples loaded: {len(items)}")
        print(f"  Y values - Mean: {np.mean(ys):.4f}, Std: {np.std(ys):.4f}")
        print(f"  Y values - Min: {np.min(ys):.4f}, Max: {np.max(ys):.4f}")
        return items, np.array(ys)
    
    def setup(self, stage=None):
        # Read data
        items, ys = self._read_items(
            self.args.uaux_csv, self.args.img_col, self.args.y_col, 
            self.args.id_col, self.args.task
        )
        
        # Train/val split
        strat = ys if (self.args.task != "regression" and self.args.stratify) else None
        idx = np.arange(len(items))
        tr_idx, va_idx = train_test_split(
            idx, test_size=self.args.val_frac, random_state=self.args.seed, stratify=strat
        )
        
        self.train_items = [items[i] for i in tr_idx]
        self.val_items = [items[i] for i in va_idx]
        self.train_ys = ys[tr_idx]
        
        # Create datasets
        self.train_dataset = Dataset(self.train_items, self.transforms)
        self.val_dataset = Dataset(self.val_items, self.transforms)
        
        # Store for logging
        self.num_train = len(self.train_items)
        self.num_val = len(self.val_items)
        # Compute stats from training set only (after train/val split)
        if self.args.task == "regression":
            self.target_stats = {
                'mean': float(np.mean(self.train_ys)),  # Use training set only
                'std': float(np.std(self.train_ys)),    # Use training set only
                'min': float(np.min(ys)),               # Keep full range for reference
                'max': float(np.max(ys))                # Keep full range for reference
            }
        else:
            self.target_stats = None
    
    def train_dataloader(self):
        if self.args.task != "regression" and self.args.balanced_sampler:
            # Class-balanced sampling
            classes, counts = np.unique(self.train_ys, return_counts=True)
            class_weights = {c: 1.0 / count for c, count in zip(classes, counts)}
            sample_weights = [class_weights[self.train_items[i]["y"]] for i in range(len(self.train_items))]
            sampler = WeightedRandomSampler(sample_weights, num_samples=len(self.train_items), replacement=True)
            
            return DataLoader(
                self.train_dataset,
                batch_size=self.args.batch_size,
                sampler=sampler,
                num_workers=self.args.num_workers,
                pin_memory=True
            )
        else:
            return DataLoader(
                self.train_dataset,
                batch_size=self.args.batch_size,
                shuffle=True,
                num_workers=self.args.num_workers,
                pin_memory=True
            )
    
    def val_dataloader(self):
        return DataLoader(
            self.val_dataset,
            batch_size=self.args.batch_size,
            shuffle=False,
            num_workers=self.args.num_workers,
            pin_memory=True
        )

def parse_args():
    ap = argparse.ArgumentParser("Train 3D ResNet on U_aux with PyTorch Lightning")
    
    # Data
    ap.add_argument("--uaux_csv", required=True, help="CSV with columns: path, y[, id]")
    ap.add_argument("--img_col", default="path", help="Image path column")
    ap.add_argument("--y_col", default="y", help="Label column")
    ap.add_argument("--id_col", default=None, help="Optional ID column")
    ap.add_argument("--val_frac", type=float, default=0.1, help="Validation fraction")
    ap.add_argument("--stratify", action="store_true", help="Stratify split (classification only)")
    ap.add_argument("--balanced_sampler", action="store_true", help="Use class-balanced sampler")
    
    # Task
    ap.add_argument("--task", choices=["regression","binary","multiclass"], default="regression")
    ap.add_argument("--num_classes", type=int, default=2, help="For multiclass")
    
    # Model & Training
    ap.add_argument("--model_name", default="resnet18", choices=list(MODEL_FNS.keys()))
    ap.add_argument("--pretrained", action="store_true", help="Use MedicalNet pretrained weights")
    ap.add_argument("--epochs", type=int, default=30)
    ap.add_argument("--batch_size", type=int, default=6)
    ap.add_argument("--lr", type=float, default=3e-4)
    ap.add_argument("--weight_decay", type=float, default=1e-3)
    ap.add_argument("--patience", type=int, default=3, help="Early stopping patience")
    ap.add_argument("--grad_clip", type=float, default=0.0, help="Gradient clipping (0 to disable)")
    ap.add_argument("--amp", action="store_true", help="Mixed precision training")
    ap.add_argument("--num_workers", type=int, default=8)
    ap.add_argument("--seed", type=int, default=1337)
    
    # Fine-tuning options
    ap.add_argument("--freeze_backbone", action="store_true",
                    help="Freeze all backbone layers (conv1..layer4); train head only unless unfreeze flags used.")
    ap.add_argument("--trainable_layers", type=str, default="fc",
                    help="Comma-separated list of layers to unfreeze: e.g. 'fc' or 'layer4,fc' or 'layer3,layer4,fc'.")
    ap.add_argument("--unfreeze_from", type=str, default=None,
                    help="Alternative to trainable_layers: unfreeze this block and all later blocks (one of: conv1,layer1,layer2,layer3,layer4,fc).")
    ap.add_argument("--lr_head", type=float, default=None,
                    help="LR for head (fc). Defaults to --lr if not set.")
    ap.add_argument("--lr_backbone", type=float, default=None,
                    help="LR for unfrozen backbone blocks. Defaults to 0.1 * --lr if not set.")
    ap.add_argument("--unfreeze_after_epochs", type=int, default=0,
                    help="Gradual unfreezing: after this many epochs, also unfreeze the next earlier block.")
    ap.add_argument("--target_norm", choices=["none", "zscore"], default="zscore",
                    help="Target normalization: 'zscore' for z-score normalization, 'none' for raw targets")
    ap.add_argument("--head_dropout", type=float, default=0.2,
                    help="Dropout rate for head layers (default: 0.2, lower than original 0.5)")
    ap.add_argument("--simple_head", action="store_true",
                    help="Use simple Linear(512->1) head instead of 2-layer MLP")
    
    # Preprocessing
    ap.add_argument("--pixdim", type=float, nargs=3, default=[1.0,1.0,1.0], help="Target spacing")
    ap.add_argument("--roi", type=int, nargs=3, default=[128,128,128], help="Target size D,H,W")
    ap.add_argument("--no_reorient", action="store_true")
    ap.add_argument("--no_resample", action="store_true")
    ap.add_argument("--no_cropfg", action="store_true")
    
    # Output
    ap.add_argument("--out_dir", required=True, help="Output directory")
    
    # Logging
    ap.add_argument("--project_name", default="brain-mri-training", help="Project name for logging")
    ap.add_argument("--use_wandb", action="store_true", help="Use Weights & Biases logging")
    ap.add_argument("--use_tensorboard", action="store_true", help="Use TensorBoard logging")
    
    # Export embeddings and predictions
    ap.add_argument("--export_csv", default=None, help="CSV to embed with trained model")
    ap.add_argument("--export_features", default=None, help="Output .npy for embeddings")
    ap.add_argument("--export_predictions", default=None, help="Output .npy for predictions")
    ap.add_argument("--export_index", default=None, help="Output .csv for embedding index")
    
    # Train/test split mode
    ap.add_argument("--use_train_test_split", action="store_true", 
                   help="Use train_labeled.csv for training and extract embeddings from test_labeled.csv")
    ap.add_argument("--test_csv", default=None, help="Path to test CSV for embedding extraction")
    
    return ap.parse_args()

def setup_loggers(args):
    """Setup PyTorch Lightning loggers"""
    loggers = []
    
    # Always use CSV logger
    csv_logger = CSVLogger(args.out_dir, name="training_logs")
    loggers.append(csv_logger)
    
    # Optional W&B logger
    if args.use_wandb and HAS_WANDB:
        wandb_logger = WandbLogger(
            project=args.project_name,
            save_dir=args.out_dir,
            name=f"{args.model_name}_{args.task}_{int(time.time())}"
        )
        loggers.append(wandb_logger)
        print("🚀 Weights & Biases logging enabled")
    
    # Optional TensorBoard logger
    if args.use_tensorboard and HAS_TENSORBOARD:
        tb_logger = TensorBoardLogger(args.out_dir, name="tensorboard_logs")
        loggers.append(tb_logger)
        print(f"📊 TensorBoard logging enabled: tensorboard --logdir {args.out_dir}/tensorboard_logs")
    
    return loggers

def main():
    args = parse_args()
    Path(args.out_dir).mkdir(parents=True, exist_ok=True)
    
    # Set seeds
    pl.seed_everything(args.seed, workers=True)
    
    print(f"🚀 Starting training with PyTorch Lightning")
    print(f"🧠 Model: {args.model_name}")
    print(f"🎯 Task: {args.task}")
    print(f"🏥 Pretrained: {'MedicalNet' if args.pretrained else 'Random init'}")
    print(f"📁 Output: {args.out_dir}")
    print(f"⚡ AMP: {args.amp}")
    
    # Handle train/test split mode
    if args.use_train_test_split:
        # Auto-configure paths for train/test split
        csv_dir = Path(args.uaux_csv).parent
        train_csv = csv_dir / "train_labeled.csv"
        test_csv = csv_dir / "test_labeled.csv"
        
        if not train_csv.exists():
            raise FileNotFoundError(f"Train CSV not found: {train_csv}. Run split_data.py first.")
        if not test_csv.exists():
            raise FileNotFoundError(f"Test CSV not found: {test_csv}. Run split_data.py first.")
            
        # Override CSV path to use train split
        args.uaux_csv = str(train_csv)
        args.test_csv = str(test_csv)
        
        # Auto-configure embedding and prediction export paths
        if not args.export_features:
            args.export_features = str(Path(args.out_dir) / "test_embeddings.npy")
        if not args.export_predictions:
            args.export_predictions = str(Path(args.out_dir) / "test_predictions.npy")
        if not args.export_index:
            args.export_index = str(Path(args.out_dir) / "test_embeddings_index.csv")
            
        print(f"🔀 Train/test split mode enabled:")
        print(f"   Training on: {train_csv}")
        print(f"   Test embeddings from: {test_csv}")
        print(f"   Embeddings output: {args.export_features}")
        print(f"   Predictions output: {args.export_predictions}")
    
    # Setup data module
    data_module = BrainMRIDataModule(args)
    data_module.setup()
    
    print(f"📊 Data: {data_module.num_train} train, {data_module.num_val} val")
    if data_module.target_stats:
        stats = data_module.target_stats
        print(f"📈 Target stats: μ={stats['mean']:.2f}, σ={stats['std']:.2f}, range=[{stats['min']:.2f}, {stats['max']:.2f}]")
        # Debug: Check if bias initialization will work
        if abs(stats['mean']) > 10:
            print(f"⚠️  WARNING: Target mean is {stats['mean']:.2f}, which is far from 0!")
            print(f"   This suggests the target statistics computation might be wrong.")
    
    # Setup model
    model = BrainMRIModule(args)
    
    # Pass target stats to model for z-scoring and bias initialization
    if data_module.target_stats and args.task == "regression":
        model.register_buffer("y_mean_buf", torch.tensor([data_module.target_stats["mean"]], dtype=torch.float32))
        model.register_buffer("y_std_buf", torch.tensor([data_module.target_stats["std"]], dtype=torch.float32))
        
        # Initialize final bias to training mean
        if model._init_bias:
            model._init_head_bias()
    total_params = sum(p.numel() for p in model.parameters())
    trainable_params = sum(p.numel() for p in model.parameters() if p.requires_grad)
    print(f"🔧 Parameters: {total_params:,} total, {trainable_params:,} trainable")
    
    # Setup callbacks
    callbacks = []
    
    # Early stopping
    early_stop = EarlyStopping(
        monitor='val_loss_epoch',
        patience=args.patience,
        mode='min',
        verbose=True
    )
    callbacks.append(early_stop)
    
    # Model checkpointing
    checkpoint = ModelCheckpoint(
        dirpath=args.out_dir,
        filename='best_model',
        monitor='val_loss_epoch',
        mode='min',
        save_top_k=1,
        verbose=True
    )
    callbacks.append(checkpoint)
    
    # Rich progress bar
    progress_bar = RichProgressBar()
    callbacks.append(progress_bar)
    
    # Setup loggers
    loggers = setup_loggers(args)
    
    # Setup trainer
    trainer = pl.Trainer(
        max_epochs=args.epochs,
        accelerator='gpu' if torch.cuda.is_available() else 'cpu',
        devices=1,
        precision='16-mixed' if args.amp else 32,
        gradient_clip_val=args.grad_clip if args.grad_clip > 0 else None,
        callbacks=callbacks,
        logger=loggers,
        deterministic='warn',  # Use 'warn' instead of True to avoid CUDA determinism issues
        enable_progress_bar=True,
        log_every_n_steps=10,
    )
    
    # Train model
    print("\n🏃 Starting training...")
    start_time = time.time()
    trainer.fit(model, data_module)
    training_time = time.time() - start_time
    
    print(f"\n🎉 Training completed in {training_time/60:.1f} minutes")
    print(f"🏆 Best validation loss: {checkpoint.best_model_score:.4f}")
    
    # Save final config
    config = vars(args)
    config.update({
        'best_val_loss': float(checkpoint.best_model_score),
        'training_time_minutes': training_time / 60,
        'total_params': total_params,
        'trainable_params': trainable_params,
        'final_train_losses': model.train_losses,
        'final_val_losses': model.val_losses,
    })
    
    with open(os.path.join(args.out_dir, "config.json"), "w") as f:
        json.dump(config, f, indent=2)
    
    # Export embeddings and predictions if requested or in train/test split mode
    export_csv_path = args.export_csv or (args.test_csv if args.use_train_test_split else None)
    
    if export_csv_path and args.export_features and args.export_index:
        print("\n🔄 Exporting embeddings and predictions...")
        
        # Load best model - need to recreate with same setup including buffers
        best_model = BrainMRIModule(args, init_bias=False)  # Don't initialize bias when loading checkpoint
        
        # Load the checkpoint state dict first
        checkpoint_state = torch.load(checkpoint.best_model_path, weights_only=False)
        state_dict = checkpoint_state["state_dict"]
        
        # Register buffers before loading if they exist in the checkpoint or need to be created
        if data_module.target_stats and args.task == "regression":
            if "y_mean_buf" in state_dict and "y_std_buf" in state_dict:
                # Buffers are in checkpoint - register them so load_state_dict can populate them
                best_model.register_buffer("y_mean_buf", torch.tensor([0.0], dtype=torch.float32))  # placeholder
                best_model.register_buffer("y_std_buf", torch.tensor([1.0], dtype=torch.float32))   # placeholder
            else:
                # Buffers not in checkpoint - create them from data stats
                best_model.register_buffer("y_mean_buf", torch.tensor([data_module.target_stats["mean"]], dtype=torch.float32))
                best_model.register_buffer("y_std_buf", torch.tensor([data_module.target_stats["std"]], dtype=torch.float32))
        
        # --- after ---
        # Load strictly so the head must match; print a diagnostic if you prefer to probe first
        k = best_model.load_state_dict(state_dict, strict=False)
        print(f"[ckpt] missing: {k.missing_keys}")
        print(f"[ckpt] unexpected: {k.unexpected_keys}")
        # guard: the head must load
        assert not any(s.startswith("model.fc") or s.startswith("fc") for s in k.missing_keys), \
            "Checkpoint did not load the head; check head architecture or state_dict keys."
        # (Optionally switch to strict=True once you’ve verified)
        # best_model.load_state_dict(state_dict, strict=True)

        # Always show the z-stats used for de-normalization
        if hasattr(best_model, "y_mean_buf") and hasattr(best_model, "y_std_buf"):
            print(f"[export] y_mean={float(best_model.y_mean_buf):.3f}  y_std={float(best_model.y_std_buf):.3f}")
        
        # Move model to appropriate device
        device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        best_model = best_model.to(device)
        
        # Create dataset for export
        export_items, _ = data_module._read_items(
            export_csv_path, args.img_col, args.y_col, args.id_col, args.task
        )
        export_dataset = Dataset(export_items, data_module.transforms)
        export_loader = DataLoader(
            export_dataset, 
            batch_size=args.batch_size, 
            shuffle=False,
            num_workers=0,  # Use 0 workers for export to avoid deadlocks
            pin_memory=False  # Disable pin_memory for export
        )
        
        # Export both embeddings and predictions if predictions path is provided
        if args.export_predictions:
            best_model.export_embeddings_and_predictions(
                export_loader, args.export_features, args.export_predictions, args.export_index
            )
        else:
            best_model.export_embeddings(export_loader, args.export_features, args.export_index)
        
        if args.use_train_test_split:
            print(f"\n🎯 Test set embeddings and predictions extracted:")
            print(f"   Features: {args.export_features}")
            if args.export_predictions:
                print(f"   Predictions: {args.export_predictions}")
            print(f"   Index: {args.export_index}")
            print(f"   Test samples: {len(export_items)}")
    
    print("\n✅ All done!")
    if args.use_wandb and HAS_WANDB:
        print("🌐 View results at: https://wandb.ai")
    if args.use_tensorboard and HAS_TENSORBOARD:
        print(f"📊 View TensorBoard: tensorboard --logdir {args.out_dir}/tensorboard_logs")

if __name__ == "__main__":
    main()
