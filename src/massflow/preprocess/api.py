from __future__ import annotations

from typing import Any, Callable, Literal, Optional, Protocol, Self, Sequence

import numpy as np

from massflow.preprocess.batch_pre_fun import BatchPreprocess
from massflow.r_preprocess.adapter import CardinalAdapter

TaskScope = Literal["batch", "dataset"]
Backend = Literal["python", "cardinal"]

class _TaskRegistrar(Protocol):
    def _register_task(
        self,
        name: str,
        *,
        scope: TaskScope = "batch",
        apply_fn: Callable[..., Sequence[Any]] | Callable[..., Any],
        **kwargs: Any,
    ) -> Self:
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
        numba_max_threads: Optional[int] = None,
    ) -> Self:
        """Register baseline correction batch task."""
        return self._register_task(
            "baseline_correction",
            scope="batch",
            apply_fn=BatchPreprocess.baseline_correction_batch,
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
            numba_max_threads=numba_max_threads,
        )

    def noise_reduction(
        self,
        *,
        method: str = "ma",
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
        numba_max_threads: Optional[int] = 10,
    ) -> Self:
        """Register noise reduction batch task."""
        return self._register_task(
            "noise_reduction",
            scope="batch",
            apply_fn=BatchPreprocess.noise_reduction_batch,
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
            numba_max_threads=numba_max_threads,
        )

    def normalization(
        self,
        *,
        scale_method: str = "none",
        method: str = "tic",
        scale: float = 1.0,
        numba_max_threads: Optional[int] = None,
    ) -> Self:
        """Register normalization batch task."""
        return self._register_task(
            "normalization",
            scope="batch",
            apply_fn=BatchPreprocess.normalization_batch,
            scale_method=scale_method,
            method=method,
            scale=scale,
            numba_max_threads=numba_max_threads,
        )

    def peak_align(
        self,
        *,
        reference: np.ndarray | None = None,
        tolerance: float | None = None,
        units: str = "ppm",
        backend: Backend = "python",
        binfun: str = "median",
        binratio: float = 2.0,
        clear_memory: bool = False,
    ) -> Self:
        """Register peak alignment task."""
        if backend == "cardinal":
            return self._register_task(
                "peak_align",
                scope="dataset",
                apply_fn=CardinalAdapter.peak_align,
                reference=reference,
                tolerance=tolerance,
                units=units,
                binfun=binfun,
                binratio=binratio,
            )

        return self._register_task(
            "peak_align",
            scope="dataset",
            apply_fn=BatchPreprocess.peak_align_batch,
            reference=reference,
            tolerance=tolerance,
            units=units,
            binfun=binfun,
            binratio=binratio,
            clear_memory=clear_memory,
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
        backend: Backend = "python",
    ) -> Self:
        """Register peak picking task."""
        if backend == "cardinal":
            method = "diff" if method == "origin" else method
            return self._register_task(
                "peak_pick",
                scope="dataset",
                apply_fn=CardinalAdapter.peak_pick,
                method=method,
                snr=snr,
                return_type=return_type,
            )

        return self._register_task(
                "peak_pick",
                scope="batch",
                apply_fn=BatchPreprocess.peak_pick_batch,
                width=width,
                method=method,
                relheight=relheight,
                snr=snr,
                return_type=return_type,
                use_numba=use_numba,
            )
