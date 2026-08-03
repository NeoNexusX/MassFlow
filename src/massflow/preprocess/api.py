from __future__ import annotations

from typing import Any, Callable, Literal, Protocol, Self
import numpy as np

from massflow.preprocess.helper.peak_align_helper import reference_computer
from massflow.preprocess.flat_pre_fun import FlatPreprocess
from massflow.preprocess.pipeline.types import TaskScope
from massflow.r_preprocess.adapter import CardinalAdapter
from massflow.preprocess.helper.stat_filter_helper import filter_peaks_by_stats, IntensityMethod
from massflow.preprocess.helper.adduct_filter_helper import filter_peaks_by_adducts
from massflow.tools.logger import get_logger
logger = get_logger("massflow.preprocess.api")

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

    batch_size: int
    temp_dir: str

    def baseline_correction(
            self,
            *,
            method: str = "locmin_numba",
            smooth: str = "none",
            span: float = 0.1,
            upper: bool = False,
            width: int = 5,
            niter: int = 15,
            baseline_scale: float = 1.0,
            m: int | None = None,
            decreasing: bool = True,
            backend: Backend = "python",
        ) -> Self:
        if backend == "cardinal":
            if niter != 15 or baseline_scale != 1.0 or m is not None:
                logger.warning(
                    "niter/baseline_scale/m are only supported by the python "
                    "baseline_correction backend"
                )
            return self._register_task(
                "baseline_correction",
                scope="dataset",
                apply_fn=CardinalAdapter.baseline_correction,
                method=method,
                smooth=smooth,
                span=span,
                upper=upper,
                width=width,
                decreasing=decreasing,
                temp_dir=self.temp_dir,
            )

        return self._register_task(
            "baseline_correction",
            scope="batch",
            apply_fn=FlatPreprocess.baseline_reduction_flat,
            method=method,
            smooth=smooth,
            span=span,
            upper=upper,
            width=width,
            niter=niter,
            baseline_scale=baseline_scale,
            m=m,
            decreasing=decreasing,
        )

    def noise_reduction(
        self,
        *,
        method: str = "gaussian_numba",
        window: int = 5,
        sd: float | None = None,
        polyorder: int = 3,
        deriv: int = 0,
        delta: float = 1.0,
        backend: Backend = "python",
    ) -> Self:
        """Register a noise reduction task."""
        if backend == "cardinal":
            return self._register_task(
                "noise_reduction",
                scope="dataset",
                apply_fn=CardinalAdapter.noise_reduction,
                method=method,
                window=window,
                sd=sd,
                polyorder=polyorder,
                deriv=deriv,
                delta=delta,
                temp_dir=self.temp_dir,
            )

        return self._register_task(
            "noise_reduction",
            scope="batch",
            apply_fn=FlatPreprocess.noise_reduction_flat,
            method=method,
            window=window,
            sd=sd,
            polyorder=polyorder,
            deriv=deriv,
            delta=delta,
        )

    def normalization(
        self,
        *,
        method: str = "tic_numba",
        scale: float | None = None,
        ref: float | None = None,
        ref_tolerance: float = 0.1,
        backend: Backend = "python",
    ) -> Self:
        """Register a normalization task."""
        if backend == "cardinal":
            if ref_tolerance != 0.1:
                logger.warning(
                    "ref_tolerance is only supported by the python "
                    "normalization backend"
                )
            return self._register_task(
                "normalization",
                scope="dataset",
                apply_fn=CardinalAdapter.normalization,
                method=method,
                scale=scale,
                ref=ref,
                temp_dir=self.temp_dir,
            )

        return self._register_task(
            "normalization",
            scope="batch",
            apply_fn=FlatPreprocess.normalization_flat,
            method=method,
            scale=scale,
            ref=ref,
            ref_tolerance=ref_tolerance,
        )

    def peak_pick(
        self,
        *,
        width: int = 5,
        method: str = "mad",
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
                temp_dir=self.temp_dir,
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

    def filter_peaks_by_isotopes(
        self,
        *,
        mz_tolerance: float = 10.0,
        units: Literal["da", "ppm"] = "ppm",
        min_charge: int = 1,
        max_charge: int = 3,
        min_isotope_count: int = 3,
        max_isotope_count: int = 10,
        use_decreasing_model: bool = True,
        start_intensity_check: int = 2,
    ) -> Self:
        """Register per-spectrum isotope filtering after peak picking."""
        return self._register_task(
            "filter_peaks_by_isotopes",
            scope="batch",
            apply_fn=FlatPreprocess.isotope_filter_flat,
            mz_tolerance=mz_tolerance,
            units=units,
            min_charge=min_charge,
            max_charge=max_charge,
            min_isotope_count=min_isotope_count,
            max_isotope_count=max_isotope_count,
            use_decreasing_model=use_decreasing_model,
            start_intensity_check=start_intensity_check,
        )

    def _reference_computer(
        self,
        *,
        binfun: str = "median",
        binratio: float = 2.0,
        tolerance: float | None = None,
        units: str = "ppm",
        min_count: int | None = None,
        min_frequency: float | None = None,
    ) -> Self:
        """Register reference computer dataset task, just for peak_align."""
        return self._register_task(
            "reference_computer",
            scope="dataset",
            apply_fn=reference_computer,
            binfun=binfun,
            binratio=binratio,
            tolerance=tolerance,
            units=units,
            min_count=min_count,
            min_frequency=min_frequency,
            batch_size=self.batch_size,
        )

    def peak_align(
        self,
        *,
        reference: np.ndarray | None = None,
        tolerance: float | None = None,
        units: str = "ppm",
        binfun: str = "min",
        binratio: float = 2.0,
        min_count: int | None = None,
        min_frequency: float | None = 0.01,
        backend: Backend = "python",
    ) -> Self:
        """Register peak alignment flat task."""
        if backend == "cardinal":
            if min_count is not None or min_frequency is not None:
                logger.warning("min_count/min_frequency are only supported by the python peak_align backend")
            return self._register_task(
                "peak_align",
                scope="dataset",
                apply_fn=CardinalAdapter.peak_align,
                reference=reference,
                tolerance=tolerance,
                units=units,
                binfun=binfun,
                binratio=binratio,
                temp_dir=self.temp_dir,
            )

        if reference is None or tolerance is None:
            self._reference_computer(
                binfun=binfun,
                binratio=binratio,
                tolerance=tolerance,
                units=units,
                min_count=min_count,
                min_frequency=min_frequency,
            )

        align_tolerance = tolerance * 1e-6 if tolerance is not None and units == "ppm" else tolerance
        return self._register_task(
            "peak_align",
            scope="batch",
            apply_fn=FlatPreprocess.peak_align_flat,
            reference=reference,
            tolerance=align_tolerance,
            units=units,
        )

    def filter_peaks_by_stats(
        self,
        *,
        min_count: int | None = None,
        min_frequency: float | None = None,
        intensity_method: IntensityMethod = "mean",
        min_intensity: float | None = None,
        intensity_quantile: float | None = 0.01,
        output_path: str | None = None,
    ) -> Self:
        """Register dataset task to filter peaks by stats like mean intensity or frequency."""
        if (
            min_count is None
            and min_frequency is None
            and min_intensity is None
            and intensity_quantile is None
        ):
            raise ValueError("At least one peak filtering threshold must be provided.")

        return self._register_task(
            "filter_peaks_by_stats",
            scope="dataset",
            apply_fn=filter_peaks_by_stats,
            min_count=min_count,
            min_frequency=min_frequency,
            intensity_method=intensity_method,
            min_intensity=min_intensity,
            intensity_quantile=intensity_quantile,
            output_path=output_path,
            temp_dir=self.temp_dir,
            batch_size=self.batch_size,
        )

    def filter_peaks_by_adducts(
        self,
        *,
        polarity: Literal["positive", "negative"] = "positive",
        ppm: float = 5.0,
        mzabs: float = 0.015,
        require_quasi: bool = True,
        min_spatial_correlation: float = 0.7,
        output_path: str | None = None,
        report_path: str | None = None,
    ) -> Self:
        """Register dataset task to filter peaks by adducts."""
        return self._register_task(
            "filter_peaks_by_adducts",
            scope="dataset",
            apply_fn=filter_peaks_by_adducts,
            polarity=polarity,
            ppm=ppm,
            mzabs=mzabs,
            require_quasi=require_quasi,
            min_spatial_correlation=min_spatial_correlation,
            batch_size=self.batch_size,
            output_path=output_path,
            temp_dir=self.temp_dir,
            report_path=report_path,
            remove=True
        )
