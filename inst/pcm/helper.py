import numpy as np


def _reshape_to_vec(x):
    x_copy = x.copy()
    if x.ndim > 1 and x.shape[1] == 1:
        x_copy = x_copy.reshape(-1)
    return x_copy


def _data_check(Y, X, Z):
    """
    Checks and (if needed) reshapes Y, X, and Z.
    For numeric (np.array) inputs, it reshapes if necessary.
    For text inputs (lists of strings), it leaves them unchanged.
    """
    # For Y: if it's an np.array and has more than one dimension, reshape to vector.
    if isinstance(Y, np.ndarray):
        if Y.ndim > 1:
            Y = Y.reshape(-1)
    # For X: if it's a numpy array, ensure it is 2D.
    if isinstance(X, np.ndarray):
        if X.ndim == 1:
            X = X[:, np.newaxis]
    # For Z: if it's a numpy array, ensure it is 2D.
    if isinstance(Z, np.ndarray):
        if Z.ndim == 1:
            Z = Z[:, np.newaxis]
    # If X or Z are lists, we assume they are already in the correct format (e.g. list of strings)
    return Y, X, Z

def _split_sample(Y, X, Z, test_split=0.5, rng=np.random.default_rng()):
    """
    Splits the data into training and test sets.
    For numeric inputs (np.array), uses numpy indexing.
    For lists (text inputs), uses list comprehensions.
    """
    # Determine the number of samples from Y.
    if isinstance(Y, np.ndarray):
        nn = Y.shape[0]
    else:
        nn = len(Y)
    
    idx_tr = rng.choice(np.arange(nn), replace=False, size=int(np.ceil(nn * (1-test_split))))
    idx_te = np.setdiff1d(np.arange(nn), idx_tr, assume_unique=True)
    
    # For Y:
    if isinstance(Y, np.ndarray):
        Ytr = Y[idx_tr]
        Yte = Y[idx_te]
    else:
        Ytr = [Y[i] for i in idx_tr]
        Yte = [Y[i] for i in idx_te]
    
    # For X:
    if isinstance(X, np.ndarray):
        Xtr = X[idx_tr, :]
        Xte = X[idx_te, :]
    else:
        Xtr = [X[i] for i in idx_tr]
        Xte = [X[i] for i in idx_te]
    
    # For Z:
    if isinstance(Z, np.ndarray):
        Ztr = Z[idx_tr, :]
        Zte = Z[idx_te, :]
    else:
        Ztr = [Z[i] for i in idx_tr]
        Zte = [Z[i] for i in idx_te]
    
    return Ytr, Xtr, Ztr, Yte, Xte, Zte
def _get_valid_args(func, args_dict):
    '''
    Return dictionary without invalid function arguments.
    Modified from https://stackoverflow.com/a/196978
    '''
    validArgs = func.__code__.co_varnames[:func.__code__.co_argcount]
    return dict((key, value) for key, value in args_dict.items()
                if key in validArgs)


def _cov_to_cor(sig):
    stds = np.sqrt(np.diag(sig))
    stds_inv = 1/stds
    cor = stds_inv * sig * stds_inv
    return cor
