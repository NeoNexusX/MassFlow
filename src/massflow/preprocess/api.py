from __future__ import annotations

from typing import Any, Callable, Protocol, Self, Sequence
import numpy as np
from massflow.preprocess.batch_pre_fun import BatchPreprocess


class _TaskRegistrar(Protocol):
    def _register_task(self, name: str, apply_fn: Callable[..., Sequence[Any]], **kwargs) -> Self:
        ...


class PreprocessorAPI(_TaskRegistrar):
    """Chainable task registration APIs for `Preprocessor`."""

    def baseline_correction(
        self,
        *,
        method: str = "asls",
        smooth: str = "none",
        span: float = 0.1,
        s: float | None = 0.0,
        upper: bool = False,
        width: int = 5,
        lam: float = 1e7,
        p: float = 0.01,
        niter: int = 15,
        baseline_scale: float = 1.0,
        m: int | None = None,
        decreasing: bool = True,
    ) -> Self:
        """Register baseline correction task.

        Parameters mirror `BatchPreprocess.baseline_correction_batch` except `batch_spectra`.
        """
        return self._register_task(
            "baseline_correction",
            BatchPreprocess.baseline_correction_batch,
            method=method,
            smooth=smooth,
            span=span,
            s=s,
            upper=upper,
            width=width,
            lam=lam,
            p=p,
            niter=niter,
            baseline_scale=baseline_scale,
            m=m,
            decreasing=decreasing,
        )

    def noise_reduction(
        self,
        *,
        method: str = "ma_numba",
        window: int = 5,
        sd: float | None = None,
        sd_intensity: float | None = None,
        p: int = 2,
        coef: np.ndarray | None = None,
        polyorder: int = 3,
        deriv: int = 0,
        delta: float = 1.0,
        wavelet: str = "db4",
        threshold_mode: str = "soft",
    ) -> Self:
        """Register noise reduction task.

        Parameters mirror `BatchPreprocess.noise_reduction_batch` except `batch_spectra`.
        """
        return self._register_task(
            "noise_reduction",
            BatchPreprocess.noise_reduction_batch,
            method=method,
            window=window,
            sd=sd,
            sd_intensity=sd_intensity,
            p=p,
            coef=coef,
            polyorder=polyorder,
            deriv=deriv,
            delta=delta,
            wavelet=wavelet,
            threshold_mode=threshold_mode,
        )

    def normalization(
        self,
        *,
        scale_method: str = "none",
        method: str = "tic",
        scale: float = 1.0,
    ) -> Self:
        """Register normalization task.

        Parameters mirror `BatchPreprocess.normalization_batch` except `batch_spectra`.
        """
        return self._register_task(
            "normalization",
            BatchPreprocess.normalization_batch,
            scale_method=scale_method,
            method=method,
            scale=scale,
        )

    def peak_align(
        self,
        *,
        ref: np.ndarray,
        tolerance: float,
        units: str = "ppm",
    ) -> Self:
        """Register peak alignment task.

        Parameters mirror `BatchPreprocess.peak_align_batch` except `batch_spectra`.
        """
        return self._register_task(
            "peak_align",
            BatchPreprocess.peak_align_batch,
            ref=ref,
            tolerance=tolerance,
            units=units,
        )

    def peak_pick(
        self,
        *,
        width: int | Sequence[int] = 2,
        method: str = "origin",
        relheight: float = 0.012,
        snr: float = 2.0,
        return_type: str = "height",
        use_numba: bool = True,
    ) -> Self:
        """Register peak picking task.

        Parameters mirror `BatchPreprocess.peak_pick_batch` except `batch_spectra`.
        """
        return self._register_task(
            "peak_pick",
            BatchPreprocess.peak_pick_batch,
            width=width,
            method=method,
            relheight=relheight,
            snr=snr,
            return_type=return_type,
            use_numba=use_numba,
        )
