from typing import Optional
import numpy as np
from logger import get_logger

logger = get_logger("ms_preprocess")

def asls_baseline(y: np.ndarray, lam: float, p: float, niter: int) -> np.ndarray:
    """Asymmetric Least Squares baseline estimation"""
    y = np.asarray(y, dtype=np.float64)
    n = y.size
    if n == 0:
        return np.array([], dtype=np.float64)
    try:
        import scipy.sparse as sp
        import scipy.sparse.linalg as spla
        D = sp.diags([1, -2, 1], [0, 1, 2], shape=(n-2, n), dtype=np.float64)
        w = np.ones(n, dtype=np.float64)
        for _ in range(int(max(1, niter))):
            W = sp.diags(w, 0, shape=(n, n), dtype=np.float64)
            Z = W + lam * (D.T @ D)
            baseline = spla.spsolve(Z, w * y)
            w = p * (y > baseline) + (1.0 - p) * (y <= baseline)
        return np.asarray(baseline, dtype=np.float64)
    except ImportError:
        D = np.zeros((n-2, n), dtype=np.float64)
        for i in range(n-2):
            D[i, i] = 1.0
            D[i, i+1] = -2.0
            D[i, i+2] = 1.0
        DT_D = D.T @ D
        w = np.ones(n, dtype=np.float64)
        baseline = y.copy()
        for _ in range(int(max(1, niter))):
            W = np.diag(w)
            Z = W + lam * DT_D
            try:
                baseline = np.linalg.solve(Z, w * y)
            except np.linalg.LinAlgError:
                baseline = np.linalg.lstsq(Z, w * y, rcond=None)[0]
            w = p * (y > baseline) + (1.0 - p) * (y <= baseline)
        return baseline

def snip_baseline(y: np.ndarray,
                  m: Optional[int] = None,
                  decreasing: bool = True,
                  epsilon: float = 1e-3) -> np.ndarray:
    """SNIP baseline estimation with adaptive early-stop."""
    y = np.asarray(y, dtype=np.float64)
    n = y.size
    if n == 0:
        return np.array([], dtype=np.float64)
    m_auto = max(1, min(50, n // 10))
    m_eval = int(max(1, min(m if m is not None else m_auto, n - 1)))
    y_work = y.copy()
    p_iter = range(m_eval, 0, -1) if decreasing else range(1, m_eval + 1)
    for pval in p_iter:
        start, end = pval, n - pval
        if end <= start:
            break
        prev_y_work_slice = y_work[start:end].copy()
        filtered_slice = 0.5 * (y_work[start - pval:end - pval] + y_work[start + pval:end + pval])
        y_work[start:end] = np.minimum(y_work[start:end], filtered_slice)
        denom = np.mean(np.abs(prev_y_work_slice)) + 1e-12
        rel_change = np.mean(np.abs(y_work[start:end] - prev_y_work_slice)) / denom
        if rel_change < epsilon:
            break
    return y_work