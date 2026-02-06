import numpy as np
import time
from functools import partial
import pytest
from massflow.module.mass_spectrum_set import MassSpectrumSet
from massflow.module.ms_data_manager_imzml import MSDataManagerImzML
from massflow.preprocess.dm_pre_fun import Preprocess
from massflow.tools.logger import get_logger
pytestmark = pytest.mark.filterwarnings("ignore:This process .* is multi-threaded, use of fork():DeprecationWarning")

logger = get_logger("test_all_normalize")
BATCH_SIZE = 512


@pytest.fixture(scope="module")
def normalize_data_manager(data_file_path="Data/other/Example_read/example.imzML") -> MSDataManagerImzML:
    mass_data = MassSpectrumSet()
    dm = MSDataManagerImzML(mass_data, filepath=data_file_path)
    dm.load_full_data_from_file()

    for _ in dm.get_batch_generator(batch_size=512):
        pass

    logger.info("data pre-load finished!")

    return dm

def run_normalization(data_manager, method="tic", scale_method="none", scale=1.0, batch_size: int = BATCH_SIZE):
    Preprocess.normalization(
        data_manager=data_manager,
        scale_method=scale_method,
        method=method,
        scale=scale,
        batch_size=batch_size,
    )

def validate_normalization_data(data_manager, method="tic", scale_method="none", scale=1.0, batch_size: int = BATCH_SIZE):
    normalized_dm = Preprocess.normalization(
        data_manager=data_manager,
        scale_method=scale_method,
        method=method,
        scale=scale,
        batch_size=batch_size,
    )

    for spectrum, spectrum_new in zip(data_manager.ms, normalized_dm.ms):
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
    def test_normalize_methods_replace_data_benchmark(self, benchmark, method, normalize_data_manager):

        logger.info(f"Testing normalization with method: {method}")
        benchmark.pedantic(
            partial(run_normalization, method=method, scale_method="none", scale=1.0, batch_size=BATCH_SIZE),
            args=(normalize_data_manager,),
            rounds=1,
            iterations=1,
            warmup_rounds=1,
        )
        logger.info(f"Finished normalization with method: {method}")

    @pytest.mark.parametrize("method", ["tic", "rms", "median"])
    def test_normalize_methods_correctness(self, method, normalize_data_manager):
        validate_normalization_data(
            normalize_data_manager,
            method=method,
            scale_method="none",
            scale=1.0,
            batch_size=BATCH_SIZE,
        )