import time
import pytest
import numpy as np
from massflow.module.mass_spectrum_set import MassSpectrumSet
from massflow.module.spectrum import Spectrum
from massflow.data_manager.ms_data_manager_imzml import MSDataManagerImzML
from massflow.preprocess.spectrum_pre_fun import SpectrumPreprocess
from massflow.tools.logger import get_logger

logger = get_logger("test_baseline")


@pytest.fixture(scope="session")
def denoised_spec() -> Spectrum:
    ms = MassSpectrumSet()
    with MSDataManagerImzML(ms, filepath="data/neg-gz4.imzML") as ms_md:
        ms_md.load_head_data()
        ms_md.inspect_data()
        sp = ms[0]
        denoised = SpectrumPreprocess.noise_reduction_spectrum(
            data=sp,
            method="gaussian",
            window=11,
        )
    return denoised


def _baseline_once(spectrum: Spectrum, method: str):
    if method == "asls":
        corrected, baseline = SpectrumPreprocess.baseline_correction_spectrum(
            data=spectrum,
            method="asls",
            lam=1e7,
            p=0.01,
            niter=15,
            baseline_scale=1.0,
        )
    elif method == "locmin":
        corrected, baseline = SpectrumPreprocess.baseline_correction_spectrum(
            data=spectrum,
            method="locmin",
            smooth="none",
            span=None, #type: ignore
            upper=False,
            width=11,
            baseline_scale=1.0,
        )
    else:
        corrected, baseline = SpectrumPreprocess.baseline_correction_spectrum(
            data=spectrum,
            method="snip",
            m=50,
            decreasing=True,
            baseline_scale=1.0,
        )
    return corrected, baseline


class TestBaseline:

    @pytest.mark.parametrize("method", ["asls", "locmin", "snip"])
    @pytest.mark.benchmark(timer=time.perf_counter)
    def test_baseline_correction_invariants(self, benchmark, method, denoised_spec):
        spec = denoised_spec
        mz0 = spec.mz_list.copy()
        int0 = spec.intensity.copy()

        def bench_callable():
            return _baseline_once(spec, method)

        result = benchmark.pedantic(
            bench_callable,
            rounds=3,
            iterations=1,
            warmup_rounds=1,
        )

        corrected, baseline = result

        x = spec.intensity
        y = corrected.intensity
        b = baseline

        assert np.array_equal(spec.mz_list, mz0)
        assert np.array_equal(spec.intensity, int0)

        assert x.shape == y.shape == b.shape
        assert np.all(np.isfinite(b))
        assert np.all(y >= 0.0)

        # expected = np.maximum(x - b, 0.0)
        # assert np.allclose(y, expected, rtol=1e-12, atol=0.0)

