from typing import Union, Optional
import numpy as np
from module.ms_module import SpectrumBaseModule
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

class MSIPreprocessor():
    """
    Abstract base class for MSI data preprocessing.
    
    This class provides a framework for implementing various preprocessing techniques
    for Mass Spectrometry Imaging (MSI) data, including peak picking, normalization,
    alignment, baseline correction, and noise reduction.
    
    Supports both traditional input (mz, msroi arrays) and MSI object input for
    better integration with the MSI framework.
    
    Attributes:
        msi_object (MSI, optional): MSI object containing data and metadata
        preprocessing_params (dict): Parameters for preprocessing operations
        processed_data (np.ndarray, optional): Processed MSI data
    """

    def __init__(self):
        """
        Initialize MSI preprocessor.
        
        Args:
            msi_object (MSI, optional): MSI object containing data and metadata
            preprocessing_params (dict, optional): Parameters for preprocessing operations
            
        Raises:
            ValueError: If neither msi_object nor (mz, msroi) are provided
        """

    @staticmethod
    def baseline_correction(
        data: Union[np.ndarray, SpectrumBaseModule],
        method: str = "asls",
        lam: float = 1e7,
        p: float = 0.01,
        niter: int = 15,
        baseline_scale: float = 0.8,
        m: Optional[int] = None,
        decreasing: bool = True,
        epsilon: float = 1e-3
    ) -> tuple[Union[np.ndarray, SpectrumBaseModule], np.ndarray]:
        """
        Remove baseline drift and background signals from mass spectra.

        Supports:
        - ASLS (asymmetric least squares): robust baseline estimation with peak preservation.
        - SNIP (statistics-sensitive non-linear iterative peak-clipping): progressive clipping with adaptive early-stop.

        Args:
            data: intensity array (np.ndarray) or SpectrumBaseModule
            method: 'asls' (default) or 'snip'
            lam: ASLS smoothness parameter (1e4-1e8, higher = smoother baseline)
            p: ASLS asymmetry parameter (0.001-0.1, lower = more peak preservation)
            niter: ASLS iterations (5-20)
            baseline_scale: scale factor (0-1) applied to the estimated baseline
            m: SNIP max half-window; None -> auto(min(50, n//10))
            decreasing: SNIP iteration order; True: p=m..1 (coarse->fine), False: p=1..m
            epsilon: SNIP adaptive early-stop threshold on relative change

        Returns:
            A tuple containing:
            - np.ndarray or SpectrumBaseModule with baseline-corrected intensity.
            - np.ndarray: The estimated baseline.
        """
        
        def apply_single(intensity: np.ndarray):
            xi = np.array(intensity, dtype=np.float64, copy=True)
            xi = np.ascontiguousarray(xi)
            # Estimate baseline via module-level functions
            if method == "asls":
                baseline = asls_baseline(xi, lam=lam, p=p, niter=niter)
            elif method == "snip":
                baseline = snip_baseline(xi, m=m, decreasing=decreasing, epsilon=epsilon)
            else:
                raise ValueError("Unsupported baseline method: use 'asls' or 'snip'")
            scale = float(np.clip(baseline_scale, 0.0, 1.0))
            scaled_baseline = scale * baseline
            corrected = xi - scaled_baseline
            corrected = np.maximum(corrected, 0.0)
            return corrected, scaled_baseline

        if isinstance(data, np.ndarray):
            return apply_single(data)
        elif isinstance(data, SpectrumBaseModule):
            corrected_intensity, baseline = apply_single(np.array(data.intensity, dtype=np.float64))
            corrected_spectrum = SpectrumBaseModule(
                mz_list=np.array(data.mz_list, dtype=np.float64, copy=True) if data.mz_list is not None else None,
                intensity=np.ascontiguousarray(corrected_intensity),
                coordinates=data.coordinates
            )
            return corrected_spectrum, baseline
        else:
            raise TypeError("data must be np.ndarray or SpectrumBaseModule")