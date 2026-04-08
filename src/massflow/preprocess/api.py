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

    def baseline_correction(self, **kwargs: Any) -> Self:
        return self._register_task(
            "baseline_correction",
            scope="batch",
            apply_fn=FlatPreprocess.baseline_correction_flat,
            **kwargs,
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

    def normalization(self, **kwargs: Any) -> Self:
        return self._register_task(
            "normalization",
            scope="batch",
            apply_fn=FlatPreprocess.normalization_flat,
            **kwargs,
        )

    def peak_pick(self, *, backend: Backend = "python", **kwargs: Any) -> Self:
        if backend == "cardinal":
            return self._register_task(
                "peak_pick",
                scope="dataset",
                apply_fn=CardinalAdapter.peak_pick,
                **kwargs,
            )

        return self._register_task(
            "peak_pick",
            scope="batch",
            apply_fn=FlatPreprocess.peak_pick_flat,
            **kwargs,
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
            scope="batch",
            apply_fn=FlatPreprocess.peak_align_flat,
            reference=reference,
            tolerance=tolerance,
            units=units,
            binfun=binfun,
            binratio=binratio,
        )
