import numpy as np
import torch
from torch.utils.data import Dataset, random_split
from transformers import (
    DistilBertForSequenceClassification,
    DistilBertTokenizerFast,
    DistilBertConfig,
    Trainer,
    TrainingArguments, 
    EarlyStoppingCallback
)
import sys
import os
# Make sure the current directory is in the path so we can import the regression module
script_dir = os.path.dirname(os.path.abspath(__file__))
if script_dir not in sys.path:
    sys.path.append(script_dir)
from regression import RegressionMethod
from sklearn.model_selection import train_test_split  # Alternatively, you can use this.

# A PyTorch Dataset for text regression.
class TextRegressionDataset(Dataset):
    def __init__(self, texts, labels, tokenizer, max_length=128):
        self.encodings = tokenizer(texts, truncation=True, padding=True, max_length=max_length)
        # Ensure labels are floats.
        self.labels = [float(label) for label in labels]
    def __getitem__(self, idx):
        item = {key: torch.tensor(val[idx]) for key, val in self.encodings.items()}
        item["labels"] = torch.tensor(self.labels[idx], dtype=torch.float)
        return item
    def __len__(self):
        return len(self.labels)


# The Hugging Face text regression model, inheriting from RegressionMethod.
class HFBTextRegressor(RegressionMethod):
    """
    A regression method that uses a Hugging Face text regression model.
    It expects:
      - X: a list of text strings.
      - Y: a list of continuous targets (floats).
    This implementation uses DistilBERT with a regression head.
    
    Parameters:
      - model_name (str): name of the pre-trained model to use.
      - train_args (dict, optional): a dictionary of arguments for TrainingArguments.
           If not provided, default values are used (including early stopping parameters).
      - **kwargs: additional keyword arguments for the model's from_pretrained method.
    """
    def __init__(self, model_name="distilbert-base-uncased", train_args=None, **kwargs):
        # Create a configuration for regression.
        config = DistilBertConfig.from_pretrained(model_name)
        config.num_labels = 1
        config.problem_type = "regression"  # use regression loss
        # Load the model with this configuration.
        model = DistilBertForSequenceClassification.from_pretrained(model_name, config=config, **kwargs)
        super().__init__(model)
        # Load the tokenizer.
        self.tokenizer = DistilBertTokenizerFast.from_pretrained(model_name)
        self.resid_type = "vanilla"
        
        # Set up training arguments with early stopping defaults if not provided.
        if train_args is None:
            train_args = {
                "output_dir": "./results",
                "num_train_epochs": 20,
                "per_device_train_batch_size": 32,
                "per_device_eval_batch_size": 32,    # Added explicit eval batch size
                "weight_decay": 0.01,                # Added weight decay for regularization
                "logging_steps": 50,
                "eval_strategy": "steps",            # evaluate every fixed number of steps
                "eval_steps": 10,                    # evaluation every 10 steps
                "save_strategy": "steps",            # save best model at epoch end
                "save_steps": 10,                    # Save every 10 steps
                "save_total_limit": 2,               # Keep only the last 2 checkpoints
                "load_best_model_at_end": True,
                "metric_for_best_model": "loss",
                "greater_is_better": False,          # Explicitly set greater_is_better to False for loss
                "warmup_ratio": 0.1,                 # Add warmup to stabilize training
                "disable_tqdm": True,
                "no_cuda": False,                    # set to False if using GPU
                "seed": 42
            }
        self.train_args = train_args  # store the dictionary for reference
        self.training_args = TrainingArguments(**train_args)
    
    def fit(self, Y, X, X_val=None, Y_val=None, val_split=0.1):
        """
        Fine-tune the model on the provided texts (X) and continuous targets (Y).
        
        Parameters:
          - Y: list of continuous targets.
          - X: list of text strings.
          - X_val, Y_val: optional validation data; if not provided, a split of size val_split is made from (X,Y).
          - val_split: proportion of training data to use for evaluation if no eval data is provided.
        """
        # Convert Y to float and standardize it to improve training stability
        Y = np.array([float(y) for y in Y])
        mean_y = np.mean(Y)
        std_y = np.std(Y) if np.std(Y) > 0 else 1.0
        Y_standardized = (Y - mean_y) / std_y
        self.y_mean = mean_y
        self.y_std = std_y
        
        # Create the training dataset.
        #full_dataset = TextRegressionDataset(X, Y, self.tokenizer)
        full_dataset = TextRegressionDataset(X, Y_standardized, self.tokenizer)
        
        if X_val is None or Y_val is None:
            # Automatically split full_dataset into training and evaluation sets.
            total = len(full_dataset)
            eval_size = int(total * val_split)
            train_size = total - eval_size
            # Use random_split to split the dataset.
            train_dataset, eval_dataset = random_split(full_dataset, [train_size, eval_size])
        else:
            # If validation data is provided, standardize it using the same parameters
            Y_val = np.array([float(y) for y in Y_val])
            Y_val_standardized = (Y_val - mean_y) / std_y
            train_dataset = TextRegressionDataset(X, Y_standardized, self.tokenizer)
            eval_dataset = TextRegressionDataset(X_val, Y_val_standardized, self.tokenizer)
        
        # Add early stopping callback with patience
        callbacks = [EarlyStoppingCallback(early_stopping_patience=3)]
        
        trainer = Trainer(
            model=self.model,
            args=self.training_args,
            train_dataset=train_dataset,
            eval_dataset=eval_dataset,
            callbacks=callbacks
        )
        trainer.train()
        self.model_fitted = self.model
        return self

    def predict(self, X):
        """
        Return the predicted continuous values for the texts in X.
        """
        dummy_labels = [0.0] * len(X)
        dataset = TextRegressionDataset(X, dummy_labels, self.tokenizer)
        trainer = Trainer(model=self.model_fitted)
        outputs = trainer.predict(dataset)
        # Unstandardize the predictions
        preds = outputs.predictions.flatten()
        unstandardized_preds = (preds * self.y_std) + self.y_mean
        return unstandardized_preds

    def residuals(self, Y, X):
        """
        Compute residuals as the difference between the true targets and the predictions.
        """
        if self.model_fitted is None:
            raise ValueError("Model not fitted yet!")
        Y = np.array(Y, dtype=float)
        preds = self.predict(X)
        return Y - preds
