import numpy as np
import time
from functools import partial

import pytest
from massflow.module.mass_spectrum_set import MassSpectrumSet
from massflow.module.ms_data_manager_imzml import MSDataManagerImzML
from massflow.preprocess.spectrum_preprocess import SpectrumPreprocess
from massflow.logger import get_logger

logger = get_logger("test_all_baseline")

def pre_load_data_setup(data_file_path="data/neg-gz4.imzML"):
    mass_data = MassSpectrumSet()
    with MSDataManagerImzML(mass_data, filepath=data_file_path) as ms_md:
        ms_md.load_full_data_from_file()
        ms_md.inspect_data()
        for spectrum_data in mass_data:
            _ = spectrum_data.intensity
        logger.info("data pre-load finished!")
    return (mass_data,), {}

def replace_baseline_data(mass_spectrum: MassSpectrumSet, method="asls"):
    for spectrum in mass_spectrum:
        x = spectrum.intensity
        mz_list = spectrum.mz_list
        if method == "asls":
            corrected, baseline = SpectrumPreprocess.baseline_correction_spectrum(
                data=spectrum,
                method="asls",
                lam=1e7,
                p=0.01,
                niter=15,
            )
        elif method == "locmin":
            corrected, baseline = SpectrumPreprocess.baseline_correction_spectrum(
                data=spectrum,
                method="locmin",
                smooth="none",
                span=None,
                upper=False,
                width=11,
            )
        else:
            corrected, baseline = SpectrumPreprocess.baseline_correction_spectrum(
                data=spectrum,
                method="snip",
                m=50,
                decreasing=True,
                baseline_scale=1.0,
            )
        
        y = corrected.intensity
        assert x.shape == y.shape and y.shape == baseline.shape
        assert y.shape == mz_list.shape
        assert np.all(y >= 0.0)

class TestAllBaseline:
    @pytest.mark.parametrize("method", ["locmin", "snip"])
    @pytest.mark.benchmark(timer=time.perf_counter)
    def test_all_baseline_methods_replace_data_benchmark(self, benchmark, method):

        logger.info(f"Testing baseline correction with method: {method}")
        benchmark.pedantic(partial(replace_baseline_data, method=method),
                           rounds=1,
                           setup=pre_load_data_setup,
                           iterations=1,
                           warmup_rounds=0)

        logger.info(f"Finished baseline correction with method: {method}")