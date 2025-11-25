import time
from functools import partial

import pytest
from massflow.module.ms_module import MS
from massflow.module.ms_data_manager_imzml import MSDataManagerImzML
from massflow.logger import get_logger
from massflow.preprocess.ms_preprocess import MSIPreprocessor

logger = get_logger("test_noise_reduction")

# @pytest.fixture(scope="function")
def pre_load_data_setup(data_file_path="data/example.imzML"):
    #pre load data
    mass_data = MS()

    with MSDataManagerImzML(mass_data, filepath=data_file_path) as ms_md:
        ms_md.load_full_data_from_file()
        ms_md.inspect_data()

        # Pre-load data
        for spectrum_data in mass_data:
            _ = spectrum_data.intensity

        logger.info("data pre-load finished!")

    return (mass_data,), {}

def replace_spectrum(mass_spectrum: MS, method="bi_ns"):
    for i, spectrum in enumerate(mass_spectrum):
        mass_spectrum[i] = MSIPreprocessor.noise_reduction_spectrum(spectrum, method=method)


def replace_spectrum_data(mass_spectrum: MS, method="bi_ns"):
    for i, spectrum in enumerate(mass_spectrum):
        spectrum_new = MSIPreprocessor.noise_reduction_spectrum(spectrum, method=method)
        mass_spectrum[i].mz_list = spectrum_new.mz_list
        mass_spectrum[i].intensity = spectrum_new.intensity


class TestNoiseReduction:

    @pytest.mark.parametrize("method", ["ma", "gaussian"])
    @pytest.mark.benchmark(timer=time.perf_counter)
    def test_noise_reduction_methods_replace_data_benchmark(self, benchmark, method):
        logger.info(f"Testing noise reduction with method: {method}")

        benchmark.pedantic(partial(replace_spectrum_data,method=method),
                           rounds=10,
                           setup=pre_load_data_setup,
                           iterations=1,
                           warmup_rounds=1)

        logger.info(f"Finished noise reduction with method: {method}")

    @pytest.mark.parametrize("method", ["ma", "gaussian"])
    @pytest.mark.benchmark(timer=time.perf_counter)
    def test_noise_reduction_methods_replace_spectrum_benchmark(self, benchmark, method):
        logger.info(f"Testing noise reduction with method: {method}")

        benchmark.pedantic(partial(replace_spectrum,method=method),
                           rounds=10,
                           setup=pre_load_data_setup,
                           iterations=1,
                           warmup_rounds=1)

        logger.info(f"Finished noise reduction with method: {method}")
