# 📊 TensorBoard Data Access Guide

This guide shows you multiple ways to access and analyze the training data exported by PyTorch Lightning to TensorBoard.

## 🎯 **Quick Access Methods**

### 1. **TensorBoard Web Interface** (Easiest)
```bash
conda activate dncit-paper
tensorboard --logdir /path/to/your/out_dir/tensorboard_logs
```
Then open: http://localhost:6006

**Benefits:**
- Interactive plots and zooming
- Real-time updates during training
- Professional visualization
- Export plots as images

### 2. **CSV Data (Already Created)**
PyTorch Lightning automatically creates CSV logs:
```bash
# View the CSV logs
ls /path/to/your/out_dir/training_logs/
cat /path/to/your/out_dir/training_logs/version_X/metrics.csv
```

### 3. **Direct Python Access**
```python
import pandas as pd
from tensorboard.backend.event_processing.event_accumulator import EventAccumulator

# Read TensorBoard event files
def read_tb_data(log_dir):
    event_acc = EventAccumulator(log_dir)
    event_acc.Reload()
    
    # Get available metrics
    print("Available metrics:", event_acc.Tags()['scalars'])
    
    # Extract specific metric
    train_loss = event_acc.Scalars('train_loss_epoch')
    steps = [event.step for event in train_loss]
    values = [event.value for event in train_loss]
    
    return steps, values

# Usage
log_dir = "/path/to/tensorboard_logs/version_X"
steps, values = read_tb_data(log_dir)
```

## 📈 **What Data is Available**

From your training run, these metrics are logged:

### **Epoch-level metrics:**
- `train_loss_epoch`: Average training loss per epoch
- `val_loss_epoch`: Average validation loss per epoch  
- `loss_diff`: Validation loss - Training loss (overfitting indicator)

### **Step-level metrics:**
- `train_loss_step`: Loss for each training batch
- `val_loss_step`: Loss for each validation batch
- `epoch`: Current epoch number

### **Hyperparameters:**
- `hparams.yaml`: All training configuration saved automatically

## 🔧 **Working Examples**

### **Example 1: Load Data with Pandas**
```python
import pandas as pd
import matplotlib.pyplot as plt

# Read the CSV logs (easier than TensorBoard events)
csv_path = "/path/to/training_logs/version_X/metrics.csv"
df = pd.read_csv(csv_path)

print(df.head())
print(f"Available columns: {df.columns.tolist()}")

# Plot training curves
plt.figure(figsize=(10, 6))
plt.plot(df['epoch'], df['train_loss_epoch'], label='Train Loss')
plt.plot(df['epoch'], df['val_loss_epoch'], label='Val Loss')
plt.xlabel('Epoch')
plt.ylabel('Loss')
plt.legend()
plt.title('Training Progress')
plt.show()
```

### **Example 2: Quick Analysis**
```python
# Quick analysis of your training
def analyze_training(csv_path):
    df = pd.read_csv(csv_path)
    
    print("Training Summary:")
    print(f"Best train loss: {df['train_loss_epoch'].min():.4f}")
    print(f"Best val loss: {df['val_loss_epoch'].min():.4f}")
    print(f"Final train loss: {df['train_loss_epoch'].iloc[-1]:.4f}")
    print(f"Final val loss: {df['val_loss_epoch'].iloc[-1]:.4f}")
    
    # Check for overfitting
    final_diff = df['val_loss_epoch'].iloc[-1] - df['train_loss_epoch'].iloc[-1]
    print(f"Overfitting indicator: {final_diff:.4f}")
    
    return df

# Usage
df = analyze_training("/path/to/metrics.csv")
```

## 📁 **File Locations**

For your training run, data is stored in:

```
{out_dir}/
├── tensorboard_logs/
│   └── version_X/
│       ├── events.out.tfevents.* (TensorBoard binary data)
│       └── hparams.yaml (hyperparameters)
├── training_logs/
│   └── version_X/
│       └── metrics.csv (CSV format data)
├── best_model.ckpt (best model checkpoint)
└── config.json (training configuration)
```

## 🚀 **Recommendations**

### **For Quick Analysis:**
1. Use the CSV files: `training_logs/version_X/metrics.csv`
2. Load with pandas for easy manipulation

### **For Interactive Exploration:**
1. Use TensorBoard web interface
2. Great for real-time monitoring during training

### **For Custom Analysis:**
1. Read TensorBoard events with EventAccumulator
2. Full programmatic control over data

### **For R Users:**
```r
# Read the CSV data in R
df <- read.csv("/path/to/metrics.csv")
library(ggplot2)

ggplot(df, aes(x = epoch)) +
  geom_line(aes(y = train_loss_epoch, color = "Train")) +
  geom_line(aes(y = val_loss_epoch, color = "Validation")) +
  labs(title = "Training Progress", y = "Loss") +
  theme_minimal()
```

## ✅ **Your Data Summary**

From your test run:
- **Training completed**: 2 epochs
- **Best validation loss**: 1004.97
- **Data location**: `/home/RDC/simnacma/H:/simnacma/CITs/Application/ixi/t1/out/ixi_resnet18_age_lightning_test/`
- **TensorBoard**: `tensorboard --logdir /home/RDC/simnacma/H:/simnacma/CITs/Application/ixi/t1/out/ixi_resnet18_age_lightning_test/tensorboard_logs`

The CSV export from our script worked and contains all your training metrics in an easy-to-use format!
