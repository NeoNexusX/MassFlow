from __future__ import annotations

from typing import Any, Callable, Literal, Protocol, Self
import numpy as np

from massflow.preprocess.flat_pre_fun import FlatPreprocess
from massflow.r_preprocess.adapter import CardinalAdapter

TaskScope = Literal["batch", "dataset"]
Backend = Literal["python", "cardinal"]


class _TaskRegistrar(Protocol):
    def _register_task(
        self,
        name: str,
        *,
        scope: TaskScope = "batch",
        apply_fn: Callable[..., Any],
        **kwargs: Any,
    ) -> Self:
        ...


class PreprocessorAPI(_TaskRegistrar):
    """Chainable task registration APIs for `Preprocessor` (flat-first)."""

    def baseline_correction(
            self,
            *,
            method: str = "locmin_numba",
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
        return self._register_task(
            "baseline_correction",
            scope="batch",
            apply_fn=FlatPreprocess.baseline_reduction_flat,
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
        """Register noise reduction flat task."""
        return self._register_task(
            "noise_reduction",
            scope="batch",
            apply_fn=FlatPreprocess.noise_reduction_flat,
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
        method: str = "tic_numba",
        scale: float | None = None,
        mz_flat: np.ndarray | None = None,
        ref: float | None = None,
        ref_tolerance: float = 0.1,
    ) -> Self:
        """Register normalization batch task."""
        return self._register_task(
            "normalization",
            scope="batch",
            apply_fn=FlatPreprocess.normalization_flat,
            method=method,
            scale=scale,
            mz_flat=mz_flat,
            ref=ref,
            ref_tolerance=ref_tolerance,
        )

    def peak_pick(
        self,
        *,
        width: int = 5,
        method: str = "quantile",
        snr: float = 2.0,
        return_type: str = "height",
        prominence: float | None = None,
        relheight: float | None = None, # Note: Only used for python backend
        nbins: int = 1,
        overlap: float = 0.5,
        backend: Backend = "python",
    ) -> Self:
        if backend == "cardinal":
            return self._register_task(
                "peak_pick",
                scope="dataset",
                apply_fn=CardinalAdapter.peak_pick,
                width=width,
                method=method,
                snr=snr,
                return_type=return_type,
                prominence=prominence,
                relheight=relheight,
                nbins=nbins,
                overlap=overlap,
            )

        return self._register_task(
            "peak_pick",
            scope="batch",
            apply_fn=FlatPreprocess.peak_pick_flat,
            width=width,
            method=method,
            snr=snr,
            return_type=return_type,
            prominence=prominence,
            relheight=relheight,
            nbins=nbins,
            overlap=overlap,
        )

    def peak_align(
        self,
        *,
        reference: np.ndarray | None = None,
        tolerance: float | None = None,
        units: str = "ppm",
        binfun: str = "median",
        binratio: float = 2.0,
        backend: Backend = "python",
    ) -> Self:
        """Register peak alignment flat task."""
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
            apply_fn=FlatPreprocess.peak_align_flat,
            reference=reference,
            tolerance=tolerance,
            units=units,
            binfun=binfun,
            binratio=binratio,
        )
