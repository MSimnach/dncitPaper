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
from sklearn.model_selection import train_test_split

# A PyTorch Dataset for text classification.
class TextClassificationDataset(Dataset):
    def __init__(self, texts, labels, tokenizer, max_length=128):
        # Tokenize texts.
        self.encodings = tokenizer(texts, truncation=True, padding=True, max_length=max_length)
        # Ensure labels are integers.
        self.labels = [int(label) for label in labels]
    def __getitem__(self, idx):
        item = {key: torch.tensor(val[idx]) for key, val in self.encodings.items()}
        item["labels"] = torch.tensor(self.labels[idx], dtype=torch.long)
        return item
    def __len__(self):
        return len(self.labels)

# The Hugging Face text classification model.
class HFBTextClassifier:
    """
    A classification method that uses a Hugging Face text classification model.
    It expects:
      - X: a list of text strings.
      - Y: a list of class labels (integers).
    This implementation uses DistilBERT with a classification head.
    
    Parameters:
      - model_name (str): name of the pre-trained model to use.
      - num_labels (int): number of classes.
      - train_args (dict, optional): dictionary of arguments for TrainingArguments.
      - **kwargs: additional keyword arguments for the model's from_pretrained method.
    """
    def __init__(self, model_name="distilbert-base-uncased", num_labels=2, train_args=None, **kwargs):
        # Set up configuration for classification.
        config = DistilBertConfig.from_pretrained(model_name)
        config.num_labels = num_labels
        config.problem_type = "single_label_classification"
        # Load the model with this configuration.
        self.model = DistilBertForSequenceClassification.from_pretrained(model_name, config=config, **kwargs)
        # Load the tokenizer.
        self.tokenizer = DistilBertTokenizerFast.from_pretrained(model_name)
        self.resid_type = "vanilla"  # This can be used by your PCM method.
        
        # Set up training arguments with early stopping defaults if not provided.
        if train_args is None:
            train_args = {
                "output_dir": "./results",
                "num_train_epochs": 20,
                "per_device_train_batch_size": 32,
                "logging_steps": 50,
                "evaluation_strategy": "steps",
                "eval_steps": 10,
                "save_strategy": "steps",
                "save_steps": 10,
                "save_total_limit": 2,
                "load_best_model_at_end": True,
                "metric_for_best_model": "loss",
                "disable_tqdm": True,
                "no_cuda": False,
                "seed": 42
            }
        self.train_args = train_args
        self.training_args = TrainingArguments(**train_args)
        self.model_fitted = None

    def fit(self, Y, X, X_val=None, Y_val=None, val_split=0.1):
        """
        Fine-tune the classification model on the provided texts (X) and labels (Y).
        If evaluation data (X_val, Y_val) is not provided, a split of size val_split is made.
        """
        full_dataset = TextClassificationDataset(X, Y, self.tokenizer)
        if X_val is None or Y_val is None:
            total = len(full_dataset)
            eval_size = int(total * val_split)
            train_size = total - eval_size
            train_dataset, eval_dataset = random_split(full_dataset, [train_size, eval_size])
        else:
            train_dataset = TextClassificationDataset(X, Y, self.tokenizer)
            eval_dataset = TextClassificationDataset(X_val, Y_val, self.tokenizer)
        
        trainer = Trainer(
            model=self.model,
            args=self.training_args,
            train_dataset=train_dataset,
            eval_dataset=eval_dataset,
            callbacks=[EarlyStoppingCallback(early_stopping_patience=self.train_args.get("early_stopping_patience", 3))]
        )
        trainer.train()
        self.model_fitted = self.model
        return self

    def predict(self, X):
        """
        Return the predicted class labels for the texts in X.
        """
        # Dummy labels are not used.
        dummy_labels = [0] * len(X)
        dataset = TextClassificationDataset(X, dummy_labels, self.tokenizer)
        trainer = Trainer(model=self.model_fitted)
        outputs = trainer.predict(dataset)
        preds = outputs.predictions.argmax(axis=1)
        return preds
    
    def predict_proba(self, X):
        """
        Return the predicted class probabilities for the texts in X.
        """
        if self.model_fitted is None:
            raise ValueError("Model not fitted yet!")
            
        # Dummy labels are not used.
        dummy_labels = [0] * len(X)
        dataset = TextClassificationDataset(X, dummy_labels, self.tokenizer)
        trainer = Trainer(model=self.model_fitted)
        outputs = trainer.predict(dataset)
        # outputs.predictions are logits; apply softmax to obtain probabilities.
        logits = outputs.predictions
        probs = torch.softmax(torch.tensor(logits), dim=1).numpy()
        return probs

    def residuals(self, Y, X):
        """
        Compute residuals as the difference between the one-hot encoded true labels 
        and the predicted probabilities.
        
        Returns an array of residuals with shape (n_samples, num_labels).
        """
        if self.model_fitted is None:
            raise ValueError("Model not fitted yet!")
        
        dummy_labels = [0] * len(X)
        dataset = TextClassificationDataset(X, dummy_labels, self.tokenizer)
        trainer = Trainer(model=self.model_fitted)
        outputs = trainer.predict(dataset)
        # outputs.predictions are logits; apply softmax to obtain probabilities.
        logits = outputs.predictions
        probs = torch.softmax(torch.tensor(logits), dim=1).numpy()
        
        Y = np.array(Y)
        num_labels = probs.shape[1]
        one_hot = np.zeros((len(Y), num_labels))
        for i, label in enumerate(Y):
            one_hot[i, label] = 1
        residuals = one_hot - probs
        return residuals
