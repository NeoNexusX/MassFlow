from __future__ import annotations

import time
from typing import Literal

import numpy as np
import pytest

from massflow.data_manager import MSDataManagerImzML
from massflow.module import MassSpectrumSet
from massflow.preprocess.preprocessor import Preprocessor
from massflow.tools.logger import get_logger

logger = get_logger("massflow.test_peak_pick_parallel")

Backend = Literal["python", "cardinal"]

ROUND = 3
FILEPATH = r"/Users/dre/Desktop/data/40TopL,10TopL,30BottomL,20BottomR-profile/40TopL,10TopL,30BottomL,20BottomR-profile.imzML"
MZ_RANGES = [(100, 1000)]
RNG_SEED = 20260412
SPECTRUM_SAMPLE_SIZE = 64
MZ_MATCH_ATOL = 1e-6
INTENSITY_MATCH_RTOL = 1e-5
INTENSITY_MATCH_ATOL = 1e-8


def run_peak_pick(
    data_manager: MSDataManagerImzML,
    backend: Backend,
    method: str,
    width: int = 5,
    snr: float = 2.0,
    return_type: str = "height",
    prominence: float | None = None,
    relheight: float | None = None,
    nbins: int = 1,
    overlap: float = 0.5,
    batch_size: int = 256,
) -> MSDataManagerImzML:
    """Run peak picking with specified parameters and return the resulting data manager."""
    return (
        Preprocessor(data_manager, batch_size=batch_size)
        .peak_pick(
            width=width,
            method=method,
            snr=snr,
            return_type=return_type,
            prominence=prominence,
            relheight=relheight,
            nbins=nbins,
            overlap=overlap,
            backend=backend,
        )
        .start()
    )

def compare_pick_result(
    py_manager: MSDataManagerImzML,
    cardinal_manager: MSDataManagerImzML,
    mz_range: tuple[float, float],
) -> None:
    mz_start, mz_end = mz_range
    if len(py_manager.ms) != len(cardinal_manager.ms):
        pytest.fail(
            f"Spectrum count mismatch: Python={len(py_manager.ms)}, cardinal={len(cardinal_manager.ms)}"
        )

    def extract_valid_peaks(spectrum) -> tuple[np.ndarray, np.ndarray]:
        mz_arr = np.asarray(spectrum.mz_list, dtype=np.float64)
        intensity = np.asarray(spectrum.intensity, dtype=np.float64)
        in_range = (mz_arr >= mz_start) & (mz_arr <= mz_end)
        valid = in_range & np.isfinite(mz_arr) & np.isfinite(intensity) & (intensity > 0)
        return mz_arr[valid], intensity[valid]

    def sample_positions(rng: np.random.Generator, size: int, sample_size: int) -> np.ndarray:
        if size <= sample_size:
            return np.arange(size, dtype=np.int64)
        return np.sort(rng.choice(size, size=sample_size, replace=False).astype(np.int64))

    def nearest_match(sorted_mz: np.ndarray, query_mz: float) -> tuple[int, float]:
        pos = int(np.searchsorted(sorted_mz, query_mz))
        best_idx = -1
        best_diff = np.inf
        for candidate in (pos - 1, pos):
            if 0 <= candidate < sorted_mz.size:
                diff = abs(sorted_mz[candidate] - query_mz)
                if diff < best_diff:
                    best_idx = candidate
                    best_diff = diff
        return best_idx, best_diff

    def assert_peak_matches(
        source_name: str,
        target_name: str,
        spectrum_idx: int,
        source_mz: np.ndarray,
        source_intensity: np.ndarray,
        target_mz: np.ndarray,
        target_intensity: np.ndarray,
    ) -> None:
        if source_mz.size == 0:
            return

        if target_mz.size == 0:
            pytest.fail(
                f"Spectrum {spectrum_idx}: {target_name} has no peaks in {mz_range}, "
                f"but {source_name} has {source_mz.size}"
            )

        if source_mz.size != target_mz.size:
            pytest.fail(
                f"Spectrum {spectrum_idx}: peak count mismatch for sampled spectrum "
                f"({source_name}={source_mz.size}, {target_name}={target_mz.size})"
            )

        for peak_idx in range(source_mz.size):
            match_idx, mz_diff = nearest_match(target_mz, source_mz[peak_idx])
            if match_idx < 0 or mz_diff > MZ_MATCH_ATOL:
                pytest.fail(
                    f"Spectrum {spectrum_idx}: {source_name} peak mz={source_mz[peak_idx]:.12f} "
                    f"has no {target_name} match within atol={MZ_MATCH_ATOL}"
                )

            if not np.isclose(
                source_intensity[peak_idx],
                target_intensity[match_idx],
                rtol=INTENSITY_MATCH_RTOL,
                atol=INTENSITY_MATCH_ATOL,
            ):
                pytest.fail(
                    f"Spectrum {spectrum_idx}: matched peak intensity mismatch for mz={source_mz[peak_idx]:.12f}\n"
                    f"{source_name}={source_intensity[peak_idx]:.12e}, "
                    f"{target_name}={target_intensity[match_idx]:.12e}"
                )

    rng = np.random.default_rng(RNG_SEED)
    sampled_spectra = sample_positions(rng, len(py_manager.ms), SPECTRUM_SAMPLE_SIZE)

    for spectrum_idx in sampled_spectra:
        py_spectrum = py_manager.ms[int(spectrum_idx)]
        cardinal_spectrum = cardinal_manager.ms[int(spectrum_idx)]
        try:
            py_mz, py_intensity = extract_valid_peaks(py_spectrum)
            cardinal_mz, cardinal_intensity = extract_valid_peaks(cardinal_spectrum)

            assert_peak_matches(
                "python",
                "cardinal",
                int(spectrum_idx),
                py_mz,
                py_intensity,
                cardinal_mz,
                cardinal_intensity,
            )
            assert_peak_matches(
                "cardinal",
                "python",
                int(spectrum_idx),
                cardinal_mz,
                cardinal_intensity,
                py_mz,
                py_intensity,
            )

            logger.info(
                "Spectrum %d full peak-match check passed: python=%d, cardinal=%d",
                int(spectrum_idx),
                py_mz.size,
                cardinal_mz.size,
            )
        finally:
            py_spectrum.clear_data()
            cardinal_spectrum.clear_data()


class TestPeakPick:
    """Test peak picking results consistency between Python and Cardinal backends."""

    @pytest.fixture(scope="class")
    def data_manager(self) -> MSDataManagerImzML:
        """prepare a data manager with raw spectra loaded for testing."""
        mass_data = MassSpectrumSet()
        dm_raw = MSDataManagerImzML(mass_data, filepath=FILEPATH)
        dm_raw.load_head_data()
        return dm_raw

    @pytest.mark.parametrize("method", ["quantile", "sd", "diff", "mad"])
    @pytest.mark.parametrize("snr", [2.0])
    @pytest.mark.parametrize("backend", ["python"])
    @pytest.mark.benchmark(timer=time.perf_counter)
    def test_peak_pick_benchmark(
        self,
        benchmark,
        data_manager: MSDataManagerImzML, # pylint: disable=redefined-outer-name
        method: str,
        snr: float,
        backend: Backend,
        width: int = 5,
        prominence: float | None = None,
        relheight: float | None = None,
        nbins: int = 1,
        overlap: float = 0.5,
        return_type: str = "height",
        batch_size: int = 256,
    ) -> None:
        """Performance Test: Run peak picking benchmark."""

        dm_picked = benchmark.pedantic(
            run_peak_pick,
            args=(data_manager,),
            kwargs={
                "backend": backend,
                "method": method,
                "width": width,
                "snr": snr,
                "return_type": return_type,
                "prominence": prominence,
                "relheight": relheight,
                "nbins": nbins,
                "overlap": overlap,
                "batch_size": batch_size,
            },
            rounds=ROUND,
            iterations=1,
            warmup_rounds=0)

        dm_picked.close()


    @pytest.mark.parametrize("method", ["quantile", "mad", "sd", "diff"])
    @pytest.mark.parametrize("mz_range", MZ_RANGES)
    @pytest.mark.parametrize("snr", [2.0])
    def test_peak_pick_consistency(
        self,
        data_manager: MSDataManagerImzML,
        method: str,
        mz_range: tuple[float, float],
        snr: float,
    ) -> None:
        """Test that peak picking results from Python and Cardinal backends are consistent within specified mz range."""

        picked_manager_python = run_peak_pick(
            data_manager,
            backend="python",
            method=method,
            width=5,
            snr=snr,
            return_type="height",
            relheight=None,
            nbins=1,
            overlap=0.5,
        )
        picked_manager_cardinal: MSDataManagerImzML | None = None

        try:
            picked_manager_cardinal = run_peak_pick(
                data_manager,
                backend="cardinal",
                method=method,
                width=5,
                snr=snr,
                return_type="height",
                relheight=None,
                nbins=1,
                overlap=0.5,
            )

            compare_pick_result(
                py_manager=picked_manager_python,
                cardinal_manager=picked_manager_cardinal,
                mz_range=mz_range,
            )
        finally:
            picked_manager_python.close()
            if picked_manager_cardinal is not None:
                picked_manager_cardinal.close()
