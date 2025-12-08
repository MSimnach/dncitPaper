############################################################
# Unified Projected Covariance Measure Test
# Combines categorical Y handling (one-hot encoding) with 
# textual/numerical input handling for X and Z
############################################################

import numpy as np
from sklearn.linear_model import LogisticRegressionCV, ElasticNetCV
from sklearn.preprocessing import LabelEncoder
from scipy import stats
from scipy.optimize import root_scalar
import logging
from typing import Dict, Callable, Tuple, Union
import sys
import os

# Add helper directory to path
from helper import _data_check, _split_sample

# Import regression methods
try:
    from .regression import RegressionMethod, RF
    from .text_regression import HFBTextRegressor
except ImportError:
    # Fallback for direct execution
    from regression import RegressionMethod, RF
    from text_regression import HFBTextRegressor
import itertools


class MultinomRegression:
    """
    Python equivalent of R's multinom_reg function using LogisticRegressionCV.
    Handles multinomial regression with LASSO/ElasticNet regularization.
    """
    
    def __init__(self, alpha: float = 0.5, nlambda: int = 40, random_state: int = 42):
        self.alpha = alpha  # ElasticNet mixing parameter (0=Ridge, 1=Lasso)
        self.nlambda = nlambda
        self.random_state = random_state
        self.model = None
        self.label_encoder = None
        self.baseline_class = None
        
    def fit(self, X: np.ndarray, y: np.ndarray):
        """Fit multinomial regression model."""
        X = np.asarray(X)
        y = np.asarray(y)
        
        # Encode labels and store baseline
        self.label_encoder = LabelEncoder()
        y_encoded = self.label_encoder.fit_transform(y)
        self.baseline_class = 0  # First class as baseline (like R's contr.treatment)
        
        # Fit LogisticRegressionCV with LASSO (simpler approach)
        self.model = LogisticRegressionCV(
            Cs=self.nlambda,  # Number of regularization strengths
            penalty='l1',  # Use LASSO penalty
            solver='liblinear',  # Solver that supports l1 penalty
            cv=10,  # Reduce CV folds for speed
            random_state=self.random_state,
            max_iter=4000  # Increase max iterations
        )
        self.model.fit(X, y_encoded)
        return self
        
    def predict_proba(self, X: np.ndarray, drop_baseline: bool = True) -> np.ndarray:
        """
        Predict probabilities, dropping baseline class by default.
        
        Returns:
        --------
        np.ndarray: (n_samples, n_classes-1) if drop_baseline=True, else (n_samples, n_classes)
        """
        X = np.asarray(X)
        probs = self.model.predict_proba(X)
        
        if drop_baseline and probs.shape[1] > 2:
            # For multiclass, drop the baseline (first) class
            probs = probs[:, 1:]
        elif drop_baseline and probs.shape[1] == 2:
            # For binary, keep only the positive class probability
            probs = probs[:, 1:2]
            
        return probs
    
    def residuals(self, X: np.ndarray, y: np.ndarray) -> np.ndarray:
        """
        Compute residuals: Y_matrix - predicted_probabilities.
        
        Returns:
        --------
        np.ndarray: (n_samples, n_classes-1) residual matrix
        """
        X = np.asarray(X)
        y = np.asarray(y)
        
        # Get predicted probabilities (drop baseline)
        probs = self.predict_proba(X, drop_baseline=True)
        
        # Create dummy matrix for y (drop baseline class)
        y_encoded = self.label_encoder.transform(y)
        n_classes = len(self.label_encoder.classes_)
        
        if n_classes == 2:
            # Binary case: Y_matrix is just the binary indicator
            Y_matrix = (y_encoded == 1).astype(float).reshape(-1, 1)
        else:
            # Multiclass case: create dummy matrix, drop baseline
            Y_matrix = np.eye(n_classes)[y_encoded][:, 1:]  # Drop first column (baseline)
            
        return Y_matrix - probs


class MGaussianRegression:
    """
    Python equivalent of R's mgaussian_reg function using ElasticNetCV.
    Handles multivariate Gaussian regression with LASSO/ElasticNet regularization.
    """
    
    def __init__(self, alpha: float = 0.5, nlambda: int = 40, random_state: int = 42):
        self.alpha = alpha
        self.nlambda = nlambda
        self.random_state = random_state
        self.models = []
        
    def fit(self, X: np.ndarray, Y: np.ndarray):
        """Fit separate ElasticNet for each response variable."""
        X = np.asarray(X)
        Y = np.asarray(Y)
        
        if Y.ndim == 1:
            Y = Y.reshape(-1, 1)
            
        self.models = []
        for j in range(Y.shape[1]):
            model = ElasticNetCV(
                l1_ratio=self.alpha,
                n_alphas=self.nlambda,
                cv=10,  # Reduce CV folds for speed
                random_state=self.random_state,
                max_iter=4000  # Increase max iterations
            )
            model.fit(X, Y[:, j])
            self.models.append(model)
        return self
        
    def predict(self, X: np.ndarray) -> np.ndarray:
        """Predict using fitted models."""
        X = np.asarray(X)
        predictions = np.column_stack([model.predict(X) for model in self.models])
        return predictions
        
    def residuals(self, X: np.ndarray, Y: np.ndarray) -> np.ndarray:
        """Compute residuals: Y - predicted_Y."""
        return np.asarray(Y) - self.predict(X)


def _is_text_input(X):
    """Check if input is textual (list of strings or numpy array with string elements)."""
    if isinstance(X, list) and len(X) > 0:
        return isinstance(X[0], str)
    elif isinstance(X, np.ndarray):
        return X.dtype == np.object_ or (X.size > 0 and isinstance(X.flat[0], str))
    return False


def _is_categorical(Y):
    """Check if Y is categorical (non-numeric or has few unique values)."""
    if isinstance(Y, (list, np.ndarray)):
        Y_array = np.asarray(Y)
        if Y_array.dtype.kind in {'U', 'S', 'O'}:  # Unicode, byte string, object
            return True
        # Check if numeric but with few unique values (heuristic)
        if len(np.unique(Y_array)) < 10 and len(Y_array) > 20:
            return True
    return False


def _handle_probability_dimensions(probs, n_classes, drop_baseline=True):
    """
    Consistently handle probability matrix dimensions for both binary and multiclass cases.
    
    Parameters:
    -----------
    probs : np.ndarray
        Probability matrix from model
    n_classes : int
        Number of classes
    drop_baseline : bool
        Whether to drop baseline class
        
    Returns:
    --------
    np.ndarray: Properly dimensioned probability matrix
    """
    # Ensure probs is 2D
    if probs.ndim == 1:
        probs = probs.reshape(-1, 1)
    
    if drop_baseline:
        if n_classes == 2 and probs.shape[1] == 2:
            # Binary case: keep only positive class probability
            probs = probs[:, 1:2]
        elif n_classes > 2 and probs.shape[1] == n_classes:
            # Multiclass case: drop baseline (first) class
            probs = probs[:, 1:]
    
    return probs


def _prepare_inputs(X, Z, text_input=False, embedded_XZ=True):
    """Prepare X and Z inputs for both text and numeric cases."""
    if text_input:
        # Convert to lists if needed
        if isinstance(X, np.ndarray):
            X = X.tolist()
        if isinstance(Z, np.ndarray):
            Z = Z.tolist()
        
        # Combine X and Z for text inputs - make embedded_XZ functional
        if embedded_XZ:
            # Structured format with clear separation (keep existing format that works)
            XZ = [f"Rede: {z} Zusammenfassung: {x}" for x, z in zip(X, Z)]
        else:
            # Simple concatenation with separator (consistent with original pcm.py)
            XZ = [f"{x} [SEP] {z}" for x, z in zip(X, Z)]
    else:
        # Numeric inputs
        if embedded_XZ:
            # Standard column stacking for numeric data
            XZ = np.column_stack([X, Z])
        else:
            # Convert numeric to text-like format for text-based regressors
            # Handle both 1D and 2D arrays properly
            X_array = np.asarray(X)
            Z_array = np.asarray(Z)
            
            # Ensure 2D arrays
            if X_array.ndim == 1:
                X_array = X_array.reshape(-1, 1)
            if Z_array.ndim == 1:
                Z_array = Z_array.reshape(-1, 1)
            
            # Create text representation for each sample
            XZ = []
            for i in range(X_array.shape[0]):
                x_str = ' '.join(map(str, X_array[i]))
                z_str = ' '.join(map(str, Z_array[i]))
                XZ.append(f"{x_str} [SEP] {z_str}")
    
    return XZ


def unified_pcm_test(Y, X, Z,
                    reg_yonxz: RegressionMethod = None,
                    reg_yonz: RegressionMethod = None,
                    reg_yhatonz: RegressionMethod = None,
                    reg_vonxz: RegressionMethod = None,
                    reg_ronz: RegressionMethod = None,
                    embedded_XZ: bool = True,
                    estimate_variance: bool = True,
                    test_split: float = 0.5,
                    max_exp: int = 50,
                    alpha: float = 0.05,
                    rng=np.random.default_rng()):
    """
    Unified PCM test that handles both categorical and continuous Y,
    as well as both textual and numerical inputs for X and Z.
    
    Parameters:
    -----------
    Y : array-like
        Outcome variable (categorical or continuous)
    X : array-like
        Treatment variables (textual or numerical)
    Z : array-like
        Conditioning variables (textual or numerical)
    reg_yonxz : RegressionMethod, optional
        Regression method for Y ~ [X, Z]
    reg_yonz : RegressionMethod, optional
        Regression method for Y ~ Z
    reg_yhatonz : RegressionMethod, optional
        Regression method for Y_hat ~ Z
    reg_vonxz : RegressionMethod, optional
        Regression method for variance estimation
    reg_ronz : RegressionMethod, optional
        Regression method for residuals ~ Z
    embedded_XZ : bool
        Whether to embed X and Z together
    estimate_variance : bool
        Whether to estimate variance
    test_split : float
        Proportion of data to use for testing
    max_exp : int
        Maximum number of expansions for variance estimation
    alpha : float
        Significance level
    rng : numpy.random.Generator
        Random number generator
        
    Returns:
    --------
    dict: Test results with p-value, test statistic, and other info
    """
    
    # Data checking and sample splitting
    Y, X, Z = _data_check(Y, X, Z)
    
    # Determine if inputs are textual or categorical
    text_input = _is_text_input(X) or _is_text_input(Z)
    categorical_Y = _is_categorical(Y)
    
    # Set default regression methods based on input types
    if reg_yonxz is None:
        reg_yonxz = HFBTextRegressor() if text_input else RF()
    if reg_yonz is None:
        reg_yonz = HFBTextRegressor() if text_input else RF()
    if reg_yhatonz is None:
        reg_yhatonz = HFBTextRegressor() if text_input else RF()
    if reg_vonxz is None:
        reg_vonxz = HFBTextRegressor() if text_input else RF()
    if reg_ronz is None:
        reg_ronz = HFBTextRegressor() if text_input else RF()
    
    # Split sample
    Ytr, Xtr, Ztr, Yte, Xte, Zte = _split_sample(Y, X, Z, test_split=test_split, rng=rng)
    
    if categorical_Y:
        # Use categorical PCM approach
        return _categorical_pcm_test(Ytr, Xtr, Ztr, Yte, Xte, Zte,
                                   reg_yonxz, reg_yonz, reg_yhatonz, 
                                   reg_vonxz, reg_ronz, text_input, 
                                   embedded_XZ, alpha)
    else:
        # Use continuous PCM approach
        return _continuous_pcm_test(Ytr, Xtr, Ztr, Yte, Xte, Zte,
                                  reg_yonxz, reg_yonz, reg_yhatonz, 
                                  reg_vonxz, reg_ronz, text_input, 
                                  embedded_XZ, estimate_variance, max_exp, alpha)


def _categorical_pcm_test(Ytr, Xtr, Ztr, Yte, Xte, Zte, 
                         reg_yonxz, reg_yonz, reg_yhatonz, 
                         reg_vonxz, reg_ronz, text_input, 
                         embedded_XZ, alpha):
    """PCM test for categorical Y using one-hot encoding."""
    
    # Prepare training inputs
    XZtr = _prepare_inputs(Xtr, Ztr, text_input, embedded_XZ)
    XZte = _prepare_inputs(Xte, Zte, text_input, embedded_XZ)
    
    # Initialize common variables
    from sklearn.preprocessing import LabelEncoder
    label_encoder = LabelEncoder()
    Ytr_encoded = label_encoder.fit_transform(Ytr)
    Yte_encoded = label_encoder.transform(Yte)
    n_classes = len(label_encoder.classes_)
    
    if text_input:
        # For textual inputs with categorical Y, we need classification methods
        # Check if the provided methods can handle text classification
        try:
            try:
                from .text_classification import HFBTextClassifier
            except ImportError:
                from text_classification import HFBTextClassifier
            
            # Use text classification methods if available
            if not hasattr(reg_yonxz, 'predict_proba'):
                # If the provided method doesn't support classification, use text classifier
                reg_yonxz = HFBTextClassifier()
                reg_yonz = HFBTextClassifier()
            
            # Y encoding already done above
            
            # Use text classification approach
            reg_yonxz.fit(Y=Ytr, X=XZtr)
            reg_yonz.fit(Y=Ytr, X=Ztr)
            
            # Store models for hhat
            fitted_ghat = reg_yonxz
            fitted_mtilde = reg_yonz
            
            # Get probabilities
            if hasattr(reg_yonxz, 'predict_proba'):
                pgtilde = reg_yonxz.predict_proba(XZtr)
                pmtilde_tr = reg_yonz.predict_proba(Ztr)
                pmtilde_te = reg_yonz.predict_proba(Zte)
            else:
                # Fallback: convert to continuous predictions
                pgtilde = reg_yonxz.predict(XZtr).reshape(-1, 1)
                pmtilde_tr = reg_yonz.predict(Ztr).reshape(-1, 1)
                pmtilde_te = reg_yonz.predict(Zte).reshape(-1, 1)
            
            # Handle probability matrix dimensions consistently
            pgtilde = _handle_probability_dimensions(pgtilde, n_classes, drop_baseline=True)
            pmtilde_tr = _handle_probability_dimensions(pmtilde_tr, n_classes, drop_baseline=True)
            pmtilde_te = _handle_probability_dimensions(pmtilde_te, n_classes, drop_baseline=True)
            
        except ImportError:
            # Fallback: use continuous regression methods for categorical Y
            # This is not ideal but better than returning fixed p-values
            logging.warning("Text classification not available, using regression fallback for categorical Y")
            
            # Encode Y as numeric for regression
            from sklearn.preprocessing import LabelEncoder
            temp_encoder = LabelEncoder()
            Ytr_numeric = temp_encoder.fit_transform(Ytr)
            
            # Use regression methods
            reg_yonxz.fit(Y=Ytr_numeric, X=XZtr)
            reg_yonz.fit(Y=Ytr_numeric, X=Ztr)
            
            # Get predictions and convert to probabilities
            pred_full = reg_yonxz.predict(XZtr)
            pred_z = reg_yonz.predict(Ztr)
            pred_z_te = reg_yonz.predict(Zte)
            
            # Convert predictions to pseudo-probabilities using sigmoid
            from scipy.special import expit
            pgtilde = expit(pred_full).reshape(-1, 1)
            pmtilde_tr = expit(pred_z).reshape(-1, 1)
            pmtilde_te = expit(pred_z_te).reshape(-1, 1)
            
    else:
        # For numeric inputs with categorical Y
        # Use the provided regression methods with appropriate adaptations
        
        # Check if methods support classification
        if hasattr(reg_yonxz, 'predict_proba'):
            # Methods support classification directly
            reg_yonxz.fit(Y=Ytr, X=XZtr)
            reg_yonz.fit(Y=Ytr, X=Ztr)
            
            # Store models for hhat
            fitted_ghat = reg_yonxz
            fitted_mtilde = reg_yonz
            
            pgtilde = reg_yonxz.predict_proba(XZtr)
            pmtilde_tr = reg_yonz.predict_proba(Ztr)
            pmtilde_te = reg_yonz.predict_proba(Zte)
            
        else:
            # Use MultinomRegression as fallback for numeric inputs
            args_glm = {"alpha": 0.5, "nlambda": 50}
            
            # Fit ghat: Y ~ X + Z
            ghat = MultinomRegression(**args_glm)
            ghat.fit(XZtr, Ytr)
            pgtilde = ghat.predict_proba(XZtr, drop_baseline=True)
            
            # Fit mtilde: Y ~ Z  
            mtilde = MultinomRegression(**args_glm)
            mtilde.fit(Ztr, Ytr)
            pmtilde_tr = mtilde.predict_proba(Ztr, drop_baseline=True)
            pmtilde_te = mtilde.predict_proba(Zte, drop_baseline=True)
            
            # Store models for later use in hhat
            fitted_ghat = ghat
            fitted_mtilde = mtilde
            
        # Handle probability matrix dimensions consistently
        pgtilde = _handle_probability_dimensions(pgtilde, n_classes, drop_baseline=True)
        pmtilde_tr = _handle_probability_dimensions(pmtilde_tr, n_classes, drop_baseline=True)
        pmtilde_te = _handle_probability_dimensions(pmtilde_te, n_classes, drop_baseline=True)
    
    # Create Y matrix (dummy coded, drop baseline)
    if n_classes == 2:
        Ymat = (Ytr_encoded == 1).astype(float).reshape(-1, 1)
        Ymat_te = (Yte_encoded == 1).astype(float).reshape(-1, 1)
    else:
        Ymat = np.eye(n_classes)[Ytr_encoded][:, 1:]  # Drop baseline
        Ymat_te = np.eye(n_classes)[Yte_encoded][:, 1:]
    
    # Compute htilde_tr = pgtilde - pmtilde_tr
    htilde_tr = pgtilde - pmtilde_tr
    
    # Compute rho
    rY_tr = Ymat - pmtilde_tr
    rho = np.mean(np.sum(rY_tr * htilde_tr, axis=1))
    
    # Store models for hhat function scope 
    # We need to access the models that were fitted above
    # Let's create a simple fallback that uses the MultinomRegression directly
    
    # Create new MultinomRegression models for hhat if needed
    if not hasattr(reg_yonxz, 'predict_proba'):
        # Use MultinomRegression for hhat
        args_glm = {"alpha": 0.5, "nlambda": 50}
        fitted_ghat = MultinomRegression(**args_glm)
        fitted_ghat.fit(XZtr, Ytr)
        fitted_mtilde = MultinomRegression(**args_glm) 
        fitted_mtilde.fit(Ztr, Ytr)
    else:
        # Use the fitted classification models
        fitted_ghat = reg_yonxz
        fitted_mtilde = reg_yonz
    
    # Define hhat function - use the fitted models directly
    def hhat(x, z):
        xz = _prepare_inputs(x, z, text_input, embedded_XZ)
        
        # Always use the stored models for consistency
        if fitted_ghat is not None and fitted_mtilde is not None:
            pg = fitted_ghat.predict_proba(xz, drop_baseline=True)
            pm = fitted_mtilde.predict_proba(z, drop_baseline=True)
        else:
            # This should never happen if the code is working correctly
            error_msg = ("Categorical PCM hhat: No fitted models available. "
                        f"fitted_ghat is None: {fitted_ghat is None}, "
                        f"fitted_mtilde is None: {fitted_mtilde is None}. "
                        "This indicates a programming error in model storage.")
            logging.error(error_msg)
            raise RuntimeError(error_msg)
        
        # Handle dimensions consistently
        pg = _handle_probability_dimensions(pg, n_classes, drop_baseline=True)
        pm = _handle_probability_dimensions(pm, n_classes, drop_baseline=True)
            
        return np.sign(rho) * (pg - pm)
    
    # Step 2: Variance estimation (simplified for categorical case)
    def vhat(x, z):
        h = hhat(x, z)
        # Use simple variance estimate based on probabilities
        if text_input and hasattr(reg_yonxz, 'predict_proba'):
            p = reg_yonxz.predict_proba(_prepare_inputs(x, z, text_input, embedded_XZ))
        elif not text_input and hasattr(reg_yonxz, 'predict_proba'):
            p = reg_yonxz.predict_proba(_prepare_inputs(x, z, text_input, embedded_XZ))
        else:
            # Use stored models if available
            if fitted_ghat is not None:
                p = fitted_ghat.predict_proba(_prepare_inputs(x, z, text_input, embedded_XZ), drop_baseline=True)
            else:
                error_msg = ("Categorical PCM vhat: No fitted model available for variance estimation. "
                            "This indicates a programming error - fitted_ghat should be available.")
                logging.error(error_msg)
                raise RuntimeError(error_msg)
        
        # Handle dimensions consistently
        p = _handle_probability_dimensions(p, n_classes, drop_baseline=True)
        
        # Ensure variance is positive and not too small
        variance = p * (1 - p)
        variance = np.maximum(variance, 1e-6)  # Avoid division by zero in fhat
        return variance
    
    # Step 3: Compute test statistic
    def fhat(x, z):
        h = hhat(x, z)
        v = vhat(x, z)
        # Avoid division by zero
        v_safe = np.maximum(v, 1e-8)
        return h / v_safe
    
    # Compute fhats on test data
    fhats = fhat(Xte, Zte)
    
    # Fit mhatf: fhats ~ Z_te using provided regression method
    if hasattr(reg_ronz, 'fit'):
        reg_ronz.fit(Y=fhats.flatten() if fhats.ndim > 1 else fhats, X=Zte)
        mhatf_pred = reg_ronz.predict(X=Zte)
        if mhatf_pred.ndim == 1:
            mhatf_pred = mhatf_pred.reshape(-1, 1)
    else:
        # Fallback to MGaussianRegression
        mhatf = MGaussianRegression()
        mhatf.fit(Zte, fhats)
        mhatf_pred = mhatf.predict(Zte)
    
    # Compute test statistic components
    rY = Ymat_te - pmtilde_te
    rT = fhats - mhatf_pred
    
    # Ensure compatible dimensions for L computation
    if rY.shape != rT.shape:
        min_cols = min(rY.shape[1], rT.shape[1])
        rY = rY[:, :min_cols]
        rT = rT[:, :min_cols]
    
    L = np.sum(rY * rT, axis=1)
    
    # Final test statistic
    L_mean = np.mean(L)
    L_std = np.std(L, ddof=1)
    
    # Handle edge cases for test statistic
    if L_std == 0 or np.isnan(L_std) or np.isinf(L_std):
        # If standard deviation is zero, this indicates a statistical problem
        if abs(L_mean) < 1e-10:
            # All L values are essentially zero - this suggests a fundamental issue
            error_msg = (f"Categorical PCM: All test statistic components are zero. "
                        f"This indicates either: (1) perfect model fit, (2) identical predictions "
                        f"from full and reduced models, or (3) numerical issues. "
                        f"L_mean={L_mean:.2e}, L_std={L_std:.2e}, n_samples={len(L)}")
            logging.error(error_msg)
            raise ValueError(error_msg)
        else:
            # Non-zero mean but zero std - likely numerical instability
            error_msg = (f"Categorical PCM: Zero standard deviation with non-zero mean. "
                        f"This indicates numerical instability or degenerate test statistics. "
                        f"L_mean={L_mean:.6f}, L_std={L_std:.2e}, n_samples={len(L)}")
            logging.error(error_msg)
            raise ValueError(error_msg)
    else:
        stat = np.sqrt(len(L)) * L_mean / L_std
        p_value = 1 - stats.norm.cdf(stat)
    
    return {
        'p_value': p_value,
        'test_statistic': stat,
        'method': 'categorical_pcm',
        'reject': p_value < alpha,
        'Y_type': 'categorical'
    }


def _continuous_pcm_test(Ytr, Xtr, Ztr, Yte, Xte, Zte,
                        reg_yonxz, reg_yonz, reg_yhatonz, 
                        reg_vonxz, reg_ronz, text_input, 
                        embedded_XZ, estimate_variance, max_exp, alpha):
    """PCM test for continuous Y using regression methods."""
    
    # Prepare training inputs
    XZtr = _prepare_inputs(Xtr, Ztr, text_input, embedded_XZ)
    
    # Fit regression of Y on [X, Z] on the training data
    reg_yonxz.fit(Y=Ytr, X=XZtr)
    yhat = reg_yonxz.predict(X=XZtr)
    
    # Fit regression of Y_hat on Z
    reg_yhatonz.fit(Y=yhat, X=Ztr)
    rho = np.mean((Ytr - reg_yhatonz.predict(X=Ztr)) * yhat)
    
    def hhat(X_, Z_):
        XZ = _prepare_inputs(X_, Z_, text_input, embedded_XZ)
        htilde = reg_yonxz.predict(XZ) - reg_yhatonz.predict(Z_)
        return np.sign(rho) * htilde
    
    # Estimate variance if required
    if estimate_variance:
        sqr = reg_yonxz.residuals(X=XZtr, Y=Ytr)**2
        reg_vonxz.fit(Y=sqr, X=XZtr)
        
        def a(c):
            if text_input:
                den = np.array(reg_vonxz.predict(XZtr))
                den = np.maximum(den, 1e-10)  # Avoid division by zero
                return np.mean(sqr / (den + c)) - 1
            else:
                if embedded_XZ:
                    den = np.column_stack([reg_vonxz.predict(XZtr),
                                         np.repeat(0, XZtr.shape[0])])
                else:
                    den = np.column_stack([np.array(reg_vonxz.predict(XZtr)),
                                         np.repeat(0, np.array(XZtr).shape[0])])
                den_max = np.maximum(np.max(den, axis=1), 1e-10)  # Avoid division by zero
                return np.mean(sqr / (den_max + c)) - 1
        
        if a(0) < 0:
            chat = 0
        else:
            lwr, upr = 0, 10
            counter = 0
            while np.sign(a(lwr)) * np.sign(a(upr)) == 1:
                upr += 5
                counter += 1
                if counter > max_exp:
                    raise ValueError("Cannot compute variance estimate, try rerunning with `estimate_variance=False`.")
            chat = root_scalar(a, method="brentq", bracket=[lwr, upr]).root
        
        def vhat(X_, Z_):
            XZ = _prepare_inputs(X_, Z_, text_input, embedded_XZ)
            if text_input:
                vtemp = np.array(reg_vonxz.predict(XZ))
                vtemp = np.maximum(vtemp, 0)
            else:
                if embedded_XZ:
                    vtemp = np.max(np.column_stack([reg_vonxz.predict(XZ),
                                                  np.repeat(0, XZ.shape[0])]), axis=1)
                else:
                    vtemp = np.max(np.column_stack([np.array(reg_vonxz.predict(XZ)),
                                                  np.repeat(0, np.array(XZ).shape[0])]), axis=1)
            return vtemp + chat
    else:
        def vhat(X_, Z_):
            return 1
    
    # Regression on the test data
    def fhat(X_, Z_):
        return hhat(X_, Z_) / vhat(X_, Z_)
    
    fhats = fhat(Xte, Zte)
    
    # Fit final regressions
    reg_ronz.fit(Y=fhats, X=Zte)
    reg_yonz.fit(Y=Yte, X=Zte)
    
    # Final test
    rY = Yte - reg_yonz.predict(X=Zte)
    rT = fhats - reg_ronz.predict(X=Zte)
    L = rY * rT
    
    # Compute test statistic
    L_mean = np.mean(L)
    L_var = np.mean(L**2) - L_mean**2
    
    if L_var <= 0 or np.isnan(L_var):
        # Zero or negative variance indicates a statistical problem
        error_msg = (f"Continuous PCM: Invalid test statistic variance. "
                    f"This indicates either: (1) degenerate residuals, (2) perfect predictions, "
                    f"or (3) numerical issues. L_mean={L_mean:.6f}, L_var={L_var:.2e}, "
                    f"n_samples={len(Yte)}")
        logging.error(error_msg)
        raise ValueError(error_msg)
    else:
        stat = np.sqrt(len(Yte)) * L_mean / np.sqrt(L_var)
        if np.isnan(stat):
            error_msg = (f"Continuous PCM: Test statistic is NaN. "
                        f"L_mean={L_mean:.6f}, L_var={L_var:.6f}, n_samples={len(Yte)}")
            logging.error(error_msg)
            raise ValueError(error_msg)
        p_value = 1 - stats.norm.cdf(stat)
    
    return {
        'p_value': p_value,
        'test_statistic': stat,
        'method': 'continuous_pcm',
        'reject': p_value < alpha,
        'Y_type': 'continuous'
    }


class UnifiedPCMTest:
    """
    Unified PCM test wrapper class that handles both categorical and continuous Y,
    as well as both textual and numerical inputs for X and Z.
    """
    
    def __init__(self, random_state: int = 42):
        self.random_state = random_state
        
    def test(self, Y, X, Z, alpha: float = 0.05, **kwargs) -> Dict:
        """
        Test H0: X ⊥ Y | Z using unified PCM approach.
        
        Parameters:
        -----------
        Y : array-like
            Outcome variable (categorical or continuous)
        X : array-like
            Treatment variables (textual or numerical)
        Z : array-like
            Conditioning variables (textual or numerical)
        alpha : float
            Significance level
        **kwargs : dict
            Additional arguments passed to unified_pcm_test
            
        Returns:
        --------
        Dict with test results
        """
        
        # Set random seed for reproducibility
        np.random.seed(self.random_state)
        rng = np.random.default_rng(self.random_state)
        
        try:
            result = unified_pcm_test(Y, X, Z, alpha=alpha, rng=rng, **kwargs)
            
            return {
                'test_statistic': result['test_statistic'],
                'p_value': result['p_value'],
                'method': result['method'],
                'reject': result['reject'],
                'Y_type': result['Y_type']
            }
            
        except Exception as e:
            logging.warning(f"Unified PCM test failed: {e}")
            return {
                'test_statistic': np.nan,
                'p_value': np.nan,
                'method': 'unified_pcm',
                'error': str(e)
            }


def unified_pcm_test_function(Y, X, Z, random_state: int = 42, alpha: float = 0.05, **kwargs) -> Tuple[float, float, str]:
    """
    Standalone function interface for unified PCM test.
    
    Parameters:
    -----------
    Y : array-like
        Outcome variable (categorical or continuous)
    X : array-like
        Treatment variables (textual or numerical)
    Z : array-like
        Conditioning variables (textual or numerical)
    random_state : int
        Random seed for reproducibility
    alpha : float
        Significance level
    **kwargs : dict
        Additional arguments passed to unified_pcm_test
        
    Returns:
    --------
    tuple: (test_statistic, p_value, method)
    """
    
    np.random.seed(random_state)
    rng = np.random.default_rng(random_state)
    
    result = unified_pcm_test(Y, X, Z, alpha=alpha, rng=rng, **kwargs)
    
    return result['test_statistic'], result['p_value'], result['method'] 