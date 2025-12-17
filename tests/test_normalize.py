import time
import pytest

import numpy as np
from massflow.module.mass_spectrum_set import MassSpectrumSet
from massflow.module.spectrum import Spectrum
from massflow.module.ms_data_manager_imzml import MSDataManagerImzML
from massflow.preprocess.spectrum_preprocess import SpectrumPreprocess
from massflow.logger import get_logger

logger = get_logger("test_normalize")


@pytest.fixture(scope="session")
def denoised_spec() -> Spectrum:
    ms = MassSpectrumSet()
    with MSDataManagerImzML(ms, filepath="data/neg-gz4.imzML") as ms_md:
        ms_md.load_full_data_from_file()
        ms_md.inspect_data()
        sp = ms[0]
        denoised = SpectrumPreprocess.noise_reduction_spectrum(
            data=sp,
            method="gaussian",
            window=11,
        )
    return denoised


def _normalize_once(spectrum: Spectrum, method: str = "tic", scale_method: str = "none", scale: float = 1.0):
    out = SpectrumPreprocess.normalization_spectrum(spectrum, method=method, scale_method=scale_method, scale=scale)
    return out


class TestNormalization:

    @pytest.mark.parametrize("method", ["median", "rms", "tic"])
    @pytest.mark.benchmark(timer=time.perf_counter)
    def test_normalize_none_invariants(self, benchmark, method, denoised_spec):
        spec = denoised_spec
        mz0 = spec.mz_list.copy()
        int0 = spec.intensity.copy()

        def bench_callable():
            return _normalize_once(spec, method=method, scale_method="none", scale=1.0)

        normalized = benchmark.pedantic(
            bench_callable,
            rounds=10,
            iterations=1,
            warmup_rounds=1,
        )

        assert np.array_equal(spec.mz_list, mz0)
        assert np.array_equal(spec.intensity, int0)
    
        y = normalized.intensity
        x = spec.intensity
        corr = float(np.corrcoef(x, y)[0, 1])
        assert np.isfinite(corr) and np.isclose(corr, 1.0, rtol=1e-6, atol=0)
        
        pos = np.where(x > 0.0)[0]
        assert pos.size >= 2
        i1, i2 = int(pos[0]), int(pos[-1])   
        old_ratio = float(x[i1] / x[i2])
        new_ratio = float(y[i1] / y[i2])
        assert np.isclose(new_ratio, old_ratio, rtol=1e-6, atol=1e-8)

        if method == "tic":
            assert np.isclose(np.sum(y), 1.0, rtol=1e-6)
        elif method == "rms":
            assert np.isclose(np.sqrt(np.mean(y ** 2)), 1.0, rtol=1e-6)
        else:
            assert np.isclose(np.median(y), 1.0, rtol=1e-6)

        logger.info(
            f"method={method} corr={corr:.6f} sum={np.sum(y):.6f} rms={np.sqrt(np.mean(y**2)):.6f} med={np.median(y):.6f}"
        )

    @pytest.mark.benchmark(timer=time.perf_counter)
    def test_normalize_unit_scaling(self, benchmark, denoised_spec):
        spec = denoised_spec
        mz0 = spec.mz_list.copy()
        int0 = spec.intensity.copy()

        def bench_callable():
            return _normalize_once(spec, method="tic", scale_method="unit", scale=1.0)

        normalized = benchmark.pedantic(
            bench_callable,
            rounds=3,
            iterations=1,
            warmup_rounds=1,
        )

        assert np.array_equal(spec.mz_list, mz0)
        assert np.array_equal(spec.intensity, int0)
        y = normalized.intensity
        assert np.min(y) >= 0.0
        assert np.max(y) <= 1.0
        x = spec.intensity
        corr = float(np.corrcoef(x, y)[0, 1])
        assert np.isfinite(corr) and np.isclose(corr, 1.0, rtol=1e-6, atol=1e-8)

        logger.info(
            f"unit-scale corr={corr:.6f} min={np.min(y):.6f} max={np.max(y):.6f}"
        )


