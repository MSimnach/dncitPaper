#!/usr/bin/env python3
"""
Analyze the extracted ResNet embeddings and predictions.
This script loads both embeddings and exported predictions to compute MSE directly
without manual inference, making the analysis faster and more efficient.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import torch
from pathlib import Path
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.linear_model import LinearRegression, Ridge, Lasso
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
# Note: ResNet model imports removed since we now use exported predictions


def evaluate_resnet_from_predictions(predictions_path, index_path, test_csv_path):
    """Evaluate the trained ResNet model using exported predictions"""
    print(f"\n🤖 ResNet Direct Evaluation (using exported predictions):")
    
    # Load exported predictions
    print(f"   Loading predictions: {predictions_path}")
    resnet_predictions = np.load(predictions_path)
    
    # Load index to get subject ID mapping
    index_df = pd.read_csv(index_path)
    
    # Load test CSV to get actual ages
    test_df = pd.read_csv(test_csv_path)
    
    # Merge to get predictions with actual ages
    results_df = index_df.merge(test_df, left_on='subject_id', right_on='id')
    # ⚠️ CRITICAL: Sort by 'row' to align with predictions array order!
    results_df = results_df.sort_values('row').reset_index(drop=True)
    # Check if merge was successful
    if len(results_df) == 0:
        print(f"❌ Error: No matching subject IDs found!")
        print(f"   Index file has {len(index_df)} subjects with IDs like: {index_df['subject_id'].head(3).tolist()}")
        print(f"   Test CSV has {len(test_df)} subjects with IDs like: {test_df['id'].head(3).tolist()}")
        print(f"   Common issue: Index uses file path stems, CSV uses different ID format")
        print(f"   Solution: Use --id_col parameter in training pipeline to specify correct ID column")
        raise ValueError("Subject ID mismatch between index and test CSV")
    
    # Extract data
    subject_ids = results_df['subject_id'].values
    actual_ages = results_df['y'].values
    
    # Handle predictions shape (flatten if needed for regression)
    if resnet_predictions.ndim > 1 and resnet_predictions.shape[1] == 1:
        resnet_predictions = resnet_predictions.flatten()
    
    # Ensure we have the same number of predictions and labels
    if len(resnet_predictions) != len(actual_ages):
        print(f"❌ Error: Array length mismatch after merge!")
        print(f"   Predictions array: {len(resnet_predictions)} items")
        print(f"   Matched subjects: {len(actual_ages)} items")
        print(f"   Index file had: {len(index_df)} entries")
        print(f"   This suggests an issue with the prediction export or index generation")
        raise ValueError(f"Mismatch: {len(resnet_predictions)} predictions vs {len(actual_ages)} labels")
    
    print(f"   Loaded {len(resnet_predictions)} predictions")
    print(f"   Predictions shape: {resnet_predictions.shape}")
    
    # Calculate ResNet metrics
    resnet_mse = mean_squared_error(actual_ages, resnet_predictions)
    resnet_rmse = np.sqrt(resnet_mse)
    resnet_r2 = r2_score(actual_ages, resnet_predictions)
    
    # Calculate per-observation squared errors
    resnet_squared_errors = (actual_ages - resnet_predictions) ** 2
    
    print(f"   ResNet performance on test data:")
    print(f"   R²: {resnet_r2:.3f}")
    print(f"   RMSE: {resnet_rmse:.2f} years")
    print(f"   MSE: {resnet_mse:.2f}")
    print(f"   Mean squared error: {resnet_squared_errors.mean():.2f}")
    print(f"   Median squared error: {np.median(resnet_squared_errors):.2f}")
    
    # Create detailed results DataFrame
    resnet_results = pd.DataFrame({
        'subject_id': subject_ids,
        'actual_age': actual_ages,
        'resnet_predicted_age': resnet_predictions,
        'resnet_residual': actual_ages - resnet_predictions,
        'resnet_squared_error': resnet_squared_errors
    })
    
    return resnet_results, (resnet_r2, resnet_rmse, resnet_mse)


def load_embeddings_with_labels(embeddings_path, index_path, test_csv_path):
    """Load embeddings and merge with original labels"""
    # Load embeddings and index
    embeddings = np.load(embeddings_path)
    index_df = pd.read_csv(index_path)
    
    # Load original test data for labels
    test_df = pd.read_csv(test_csv_path)
    
    # Merge to get embeddings with labels
    results_df = index_df.merge(test_df, left_on='subject_id', right_on='id')
    # ⚠️ CRITICAL: Sort by 'row' to align with embeddings array order!
    results_df = results_df.sort_values('row').reset_index(drop=True)
    
    print(f"📊 Data loaded:")
    print(f"   Embeddings shape: {embeddings.shape}")
    print(f"   Samples with labels: {len(results_df)}")
    print(f"   Feature dimension: {embeddings.shape[1]}")
    
    return embeddings, results_df


def analyze_embeddings(embeddings, labels_df):
    """Perform basic analysis of the embeddings"""
    ages = labels_df['y'].values
    
    print(f"\n🔍 Embedding Analysis:")
    print(f"   Age range: [{ages.min():.1f}, {ages.max():.1f}]")
    print(f"   Age mean ± std: {ages.mean():.1f} ± {ages.std():.1f}")
    print(f"   Embedding mean ± std: {embeddings.mean():.3f} ± {embeddings.std():.3f}")
    print(f"   Embedding range: [{embeddings.min():.3f}, {embeddings.max():.3f}]")
    
    # Check for any NaN or infinite values
    nan_count = np.isnan(embeddings).sum()
    inf_count = np.isinf(embeddings).sum()
    print(f"   NaN values: {nan_count}")
    print(f"   Infinite values: {inf_count}")


def test_linear_prediction(embeddings, ages, subject_ids):
    """Test how well the embeddings can predict age with linear regression"""
    print(f"\n🧮 Linear Regression Test:")
    
    # Split embeddings for testing
    X_train, X_test, y_train, y_test, ids_train, ids_test = train_test_split(
        embeddings, ages, subject_ids, test_size=0.3, random_state=42
    )
    
    # Standardize features to improve regularization
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)
    
    # Test multiple models to address overfitting
    models = {
        'Linear': LinearRegression(),
        'Ridge (α=1.0)': Ridge(alpha=1.0),
        'Ridge (α=10.0)': Ridge(alpha=10.0),
        'Lasso (α=0.1)': Lasso(alpha=0.1)
    }
    
    results = {}
    best_model = None
    best_test_r2 = -np.inf
    
    print(f"   Training on {len(X_train)} samples, testing on {len(X_test)} samples")
    print(f"   Model comparison:")
    
    for name, model in models.items():
        # Fit model
        model.fit(X_train_scaled, y_train)
        
        # Predict
        y_pred_train = model.predict(X_train_scaled)
        y_pred_test = model.predict(X_test_scaled)
        
        # Calculate metrics
        train_r2 = r2_score(y_train, y_pred_train)
        test_r2 = r2_score(y_test, y_pred_test)
        train_rmse = np.sqrt(mean_squared_error(y_train, y_pred_train))
        test_rmse = np.sqrt(mean_squared_error(y_test, y_pred_test))
        
        results[name] = {
            'model': model,
            'train_r2': train_r2,
            'test_r2': test_r2,
            'train_rmse': train_rmse,
            'test_rmse': test_rmse,
            'y_pred_train': y_pred_train,
            'y_pred_test': y_pred_test
        }
        
        print(f"     {name:12}: Train R²={train_r2:.3f}, Test R²={test_r2:.3f}, Test RMSE={test_rmse:.1f}")
        
        # Track best model by test R²
        if test_r2 > best_test_r2:
            best_test_r2 = test_r2
            best_model = name
    
    print(f"   Best model: {best_model}")
    
    # Use best model for detailed analysis
    best_results = results[best_model]
    y_pred_train = best_results['y_pred_train'] 
    y_pred_test = best_results['y_pred_test']
    
    # Calculate per-observation squared errors
    train_squared_errors = (y_train - y_pred_train) ** 2
    test_squared_errors = (y_test - y_pred_test) ** 2
    
    # Create detailed results DataFrames
    train_results = pd.DataFrame({
        'subject_id': ids_train,
        'actual_age': y_train,
        'predicted_age': y_pred_train,
        'residual': y_train - y_pred_train,
        'squared_error': train_squared_errors,
        'split': 'train'
    })
    
    test_results = pd.DataFrame({
        'subject_id': ids_test,
        'actual_age': y_test,
        'predicted_age': y_pred_test,
        'residual': y_test - y_pred_test,
        'squared_error': test_squared_errors,
        'split': 'test'
    })
    
    # Combine results
    detailed_results = pd.concat([train_results, test_results], ignore_index=True)
    
    # Calculate overall metrics
    train_r2 = r2_score(y_train, y_pred_train)
    test_r2 = r2_score(y_test, y_pred_test)
    train_rmse = np.sqrt(mean_squared_error(y_train, y_pred_train))
    test_rmse = np.sqrt(mean_squared_error(y_test, y_pred_test))
    
    # Print final results for best model
    final_train_r2 = best_results['train_r2']
    final_test_r2 = best_results['test_r2'] 
    final_train_rmse = best_results['train_rmse']
    final_test_rmse = best_results['test_rmse']
    
    print(f"\n   Final results ({best_model}):")
    print(f"   Train R²: {final_train_r2:.3f}")
    print(f"   Test R²: {final_test_r2:.3f}")
    print(f"   Train RMSE: {final_train_rmse:.2f} years")
    print(f"   Test RMSE: {final_test_rmse:.2f} years")
    print(f"   Mean test squared error: {test_squared_errors.mean():.2f}")
    print(f"   Median test squared error: {np.median(test_squared_errors):.2f}")
    print(f"   Overfitting ratio (Train R² / Test R²): {final_train_r2 / max(final_test_r2, 0.001):.2f}")
    
    return best_results['model'], (y_test, y_pred_test), detailed_results


def dimensionality_reduction(embeddings, labels_df, n_components=2):
    """Perform PCA and t-SNE for visualization"""
    print(f"\n🎯 Dimensionality Reduction:")
    
    ages = labels_df['y'].values
    
    # PCA
    print("   Computing PCA...")
    pca = PCA(n_components=n_components)
    embeddings_pca = pca.fit_transform(embeddings)
    explained_var = pca.explained_variance_ratio_
    print(f"   PCA explained variance: {explained_var.sum():.3f} ({explained_var[0]:.3f}, {explained_var[1]:.3f})")
    
    # t-SNE (on PCA-reduced data for speed)
    print("   Computing t-SNE...")
    # Use first 50 PCA components for t-SNE to speed up computation
    pca_50 = PCA(n_components=min(50, embeddings.shape[1]))
    embeddings_pca_50 = pca_50.fit_transform(embeddings)
    
    tsne = TSNE(n_components=n_components, random_state=42, perplexity=30)
    embeddings_tsne = tsne.fit_transform(embeddings_pca_50)
    
    return embeddings_pca, embeddings_tsne, ages


def save_analysis_plots(embeddings_pca, embeddings_tsne, ages, lr_results, output_dir):
    """Save analysis plots"""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Create figure with subplots
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle('ResNet Brain MRI Embedding Analysis', fontsize=16)
    
    # PCA plot
    scatter = axes[0, 0].scatter(embeddings_pca[:, 0], embeddings_pca[:, 1], c=ages, cmap='viridis', alpha=0.7)
    axes[0, 0].set_title('PCA Visualization')
    axes[0, 0].set_xlabel('PC1')
    axes[0, 0].set_ylabel('PC2')
    plt.colorbar(scatter, ax=axes[0, 0], label='Age (years)')
    
    # t-SNE plot
    scatter2 = axes[0, 1].scatter(embeddings_tsne[:, 0], embeddings_tsne[:, 1], c=ages, cmap='viridis', alpha=0.7)
    axes[0, 1].set_title('t-SNE Visualization')
    axes[0, 1].set_xlabel('t-SNE 1')
    axes[0, 1].set_ylabel('t-SNE 2')
    plt.colorbar(scatter2, ax=axes[0, 1], label='Age (years)')
    
    # Age distribution
    axes[1, 0].hist(ages, bins=20, alpha=0.7, edgecolor='black')
    axes[1, 0].set_title('Age Distribution')
    axes[1, 0].set_xlabel('Age (years)')
    axes[1, 0].set_ylabel('Count')
    
    # Linear regression results
    if lr_results:
        y_test, y_pred_test = lr_results
        axes[1, 1].scatter(y_test, y_pred_test, alpha=0.7)
        axes[1, 1].plot([y_test.min(), y_test.max()], [y_test.min(), y_test.max()], 'r--', lw=2)
        axes[1, 1].set_title('Linear Regression: Predicted vs Actual Age')
        axes[1, 1].set_xlabel('Actual Age')
        axes[1, 1].set_ylabel('Predicted Age')
        
        # Add R² to the plot
        r2 = r2_score(y_test, y_pred_test)
        axes[1, 1].text(0.05, 0.95, f'R² = {r2:.3f}', transform=axes[1, 1].transAxes, 
                       bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    plt.tight_layout()
    
    # Save plot
    plot_path = output_dir / 'embedding_analysis.png'
    plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    print(f"   Saved analysis plot: {plot_path}")
    plt.close()


def main(results_dir=None, test_csv_path=None, output_dir=None):
    """
    Main analysis function
    
    Args:
        results_dir: Directory containing test_embeddings.npy, test_predictions.npy, and test_embeddings_index.csv
        test_csv_path: Path to test_labeled.csv file
        output_dir: Directory to save analysis results
    """
    # Default paths (can be overridden by arguments)
    if results_dir is None:
        results_dir = "/home/RDC/simnacma/H:/simnacma/CITs/Application/ixi/t1/results"
    if test_csv_path is None:
        test_csv_path = "/home/RDC/simnacma/H:/simnacma/CITs/Application/ixi/t1/csv/test_labeled.csv"
    if output_dir is None:
        output_dir = f"{results_dir}/analysis"
    
    # Construct file paths
    results_path = Path(results_dir)
    embeddings_path = results_path / "test_embeddings.npy"
    predictions_path = results_path / "test_predictions.npy"
    index_path = results_path / "test_embeddings_index.csv"
    
    print("🧠 ResNet Brain MRI Embedding Analysis")
    print("=" * 50)
    print(f"📁 Input files:")
    print(f"   Embeddings: {embeddings_path}")
    print(f"   Predictions: {predictions_path}")
    print(f"   Index: {index_path}")
    print(f"   Test CSV: {test_csv_path}")
    print(f"   Output: {output_dir}")
    
    # Check if predictions file exists
    if not predictions_path.exists():
        print(f"❌ Error: Predictions file not found: {predictions_path}")
        print("   Make sure you've run the training pipeline with the updated version that exports predictions.")
        print("   The file test_predictions.npy should be generated alongside test_embeddings.npy")
        return
    
    # Evaluate ResNet model using exported predictions
    resnet_results, resnet_metrics = evaluate_resnet_from_predictions(str(predictions_path), str(index_path), test_csv_path)
    
    # Load embeddings data
    embeddings, labels_df = load_embeddings_with_labels(str(embeddings_path), str(index_path), test_csv_path)
    ages = labels_df['y'].values
    
    # Basic analysis
    analyze_embeddings(embeddings, labels_df)
    
    # Test linear prediction on embeddings
    subject_ids = labels_df['subject_id'].values
    lr, lr_results, detailed_results = test_linear_prediction(embeddings, ages, subject_ids)
    
    # Dimensionality reduction
    embeddings_pca, embeddings_tsne, ages = dimensionality_reduction(embeddings, labels_df)
    
    # Merge ResNet results with embedding-based results
    print(f"\n🔄 Merging Results:")
    
    # Merge ResNet results with embedding results (both use all test data)
    combined_results = resnet_results.merge(
        labels_df[['subject_id', 'y']], 
        on='subject_id'
    ).rename(columns={'y': 'actual_age_check'})
    
    # Add embedding results for the test split from detailed_results
    test_embedding_results = detailed_results[detailed_results['split'] == 'test'].copy()
    test_embedding_results = test_embedding_results.rename(columns={
        'predicted_age': 'embedding_predicted_age',
        'residual': 'embedding_residual', 
        'squared_error': 'embedding_squared_error'
    })
    
    # Merge all results
    final_results = combined_results.merge(
        test_embedding_results[['subject_id', 'embedding_predicted_age', 'embedding_residual', 'embedding_squared_error']],
        on='subject_id',
        how='left'  # Keep all ResNet results, add embedding results where available
    )
    
    print(f"   Combined results: {len(final_results)} samples")
    print(f"   ResNet predictions: {len(final_results)} samples")
    print(f"   Embedding predictions: {final_results['embedding_predicted_age'].notna().sum()} samples")
    
    # Compare methods
    print(f"\n📈 Method Comparison:")
    resnet_r2, resnet_rmse, resnet_mse = resnet_metrics
    embedding_r2 = r2_score(lr_results[0], lr_results[1])
    embedding_rmse = np.sqrt(mean_squared_error(lr_results[0], lr_results[1]))
    embedding_mse = mean_squared_error(lr_results[0], lr_results[1])
    
    print(f"   ResNet (direct):     R²={resnet_r2:.3f}, RMSE={resnet_rmse:.1f}, MSE={resnet_mse:.1f}")
    print(f"   Embeddings + Ridge: R²={embedding_r2:.3f}, RMSE={embedding_rmse:.1f}, MSE={embedding_mse:.1f}")
    
    if resnet_r2 > embedding_r2:
        print(f"   🏆 ResNet direct prediction is better by {resnet_r2 - embedding_r2:.3f} R² points")
    else:
        print(f"   🏆 Embedding approach is better by {embedding_r2 - resnet_r2:.3f} R² points")
    
    # Save analysis plots
    print(f"\n📊 Saving Analysis:")
    save_analysis_plots(embeddings_pca, embeddings_tsne, ages, lr_results, output_dir)
    
    # Save detailed per-observation results (now with both ResNet and embedding predictions)
    results_csv_path = Path(output_dir) / 'per_observation_results.csv'
    final_results.to_csv(results_csv_path, index=False)
    print(f"   Saved per-observation results: {results_csv_path}")
    
    # Save embedding-only results for backward compatibility
    embedding_results_path = Path(output_dir) / 'embedding_only_results.csv'
    detailed_results.to_csv(embedding_results_path, index=False)
    print(f"   Saved embedding-only results: {embedding_results_path}")
    
    # Save processed data
    analysis_data = {
        'embeddings_pca': embeddings_pca,
        'embeddings_tsne': embeddings_tsne,
        'ages': ages,
        'subject_ids': labels_df['subject_id'].values
    }
    
    analysis_path = Path(output_dir) / 'analysis_data.npz'
    np.savez(analysis_path, **analysis_data)
    print(f"   Saved analysis data: {analysis_path}")
    
    print(f"\n✅ Analysis completed!")
    print(f"   Results saved in: {output_dir}")
    print(f"\n💡 Summary:")
    print(f"   - {len(embeddings)} test samples with 512-dimensional embeddings")
    print(f"   - ResNet direct R² = {resnet_r2:.3f} (MSE = {resnet_mse:.1f})")
    print(f"   - Embedding Ridge R² = {embedding_r2:.3f} (MSE = {embedding_mse:.1f})")
    print(f"   - Per-observation MSE saved for both ResNet and embedding methods")
    print(f"   - ResNet embeddings capture meaningful age-related features from brain MRIs")


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Analyze ResNet embeddings and predictions")
    parser.add_argument("--results_dir", default=None,
                       help="Directory containing test_embeddings.npy, test_predictions.npy, and test_embeddings_index.csv")
    parser.add_argument("--test_csv", default=None,
                       help="Path to test_labeled.csv file")
    parser.add_argument("--output_dir", default=None,
                       help="Directory to save analysis results")
    
    args = parser.parse_args()
    
    main(results_dir=args.results_dir, test_csv_path=args.test_csv, output_dir=args.output_dir)
