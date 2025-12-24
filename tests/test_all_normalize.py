import numpy as np
import time
from functools import partial
import pytest
from massflow.module.mass_spectrum_set import MassSpectrumSet
from massflow.module.ms_data_manager_imzml import MSDataManagerImzML
from massflow.preprocess.spectrum_preprocess import SpectrumPreprocess
from massflow.logger import get_logger

logger = get_logger("test_all_normalize")

def pre_load_data_setup(data_file_path="data/neg-gz4.imzML"):
    mass_data = MassSpectrumSet()

    with MSDataManagerImzML(mass_data, filepath=data_file_path) as ms_md:
        ms_md.load_full_data_from_file()
        ms_md.inspect_data()

        # Pre-load data
        for spectrum_data in mass_data:
            _ = spectrum_data.intensity

        logger.info("data pre-load finished!")

    return (mass_data,), {}


def replace_normalization_data(mass_spectrum: MassSpectrumSet, method="tic", scale_method="none", scale=1.0):
    for spectrum in mass_spectrum:
        spectrum_new = SpectrumPreprocess.normalization_spectrum(
            spectrum, method=method, scale_method=scale_method, scale=scale
        )

        x = spectrum.intensity
        y = spectrum_new.intensity

        assert x.shape == y.shape
        if method == "tic":
            assert np.isclose(np.sum(y), 1, rtol=1e-6)
        elif method == "rms":
            assert np.isclose(np.sqrt(np.mean(y ** 2)), 1.0, rtol=1e-6)
        else:
            assert np.isclose(np.median(y), 1.0, rtol=1e-6)
        if len(x) >= 2 and len(y) >= 2:
            corr = float(np.corrcoef(x, y)[0, 1])
            assert np.isfinite(corr) and corr > 0.999
            pos = np.where(x > 0.0)[0]
            if pos.size >= 2:
                i1, i2 = int(pos[0]), int(pos[-1])
                old_ratio = float(x[i1] / x[i2])
                new_ratio = float(y[i1] / y[i2])
                assert np.isclose(new_ratio, old_ratio, rtol=1e-6, atol=1e-8)


class TestAllNormalization:
    @pytest.mark.parametrize("method", ["tic", "rms", "median"])
    @pytest.mark.benchmark(timer=time.perf_counter)
    def test_normalize_methods_replace_data_benchmark(self, benchmark, method):

        logger.info(f"Testing normalization with method: {method}")
        benchmark.pedantic(
            partial(replace_normalization_data, method=method, scale_method="none", scale=1.0),
            setup=pre_load_data_setup,
            rounds=10,
            iterations=1,
            warmup_rounds=1,
        )
        logger.info(f"Finished normalization with method: {method}")