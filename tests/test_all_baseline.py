import numpy as np
import time
from functools import partial
import pytest
from massflow.module.mass_spectrum_set import MassSpectrumSet
from massflow.module.ms_data_manager_imzml import MSDataManagerImzML
from massflow.preprocess.dm_pre_fun import Preprocess
from massflow.tools.logger import get_logger
pytestmark = pytest.mark.filterwarnings("ignore:This process .* is multi-threaded, use of fork():DeprecationWarning")

logger = get_logger("test_all_baseline")
BATCH_SIZE = 512

@pytest.fixture(scope="module")
def baseline_data_manager(data_file_path="Data/other/Example_read/example.imzML") -> MSDataManagerImzML:
    mass_data = MassSpectrumSet()
    dm = MSDataManagerImzML(mass_data, filepath=data_file_path)
    dm.load_full_data_from_file()
    #dm.inspect_data()
    for _ in dm.get_batch_generator(batch_size=512):
        pass
    logger.info("data pre-load finished!")
    return dm

def run_baseline_correction(data_manager: MSDataManagerImzML, method="asls", batch_size: int = BATCH_SIZE):
    if method == "asls":
        Preprocess.baseline_correction(
            data_manager=data_manager,
            method="asls",
            lam=1e7,
            p=0.01,
            niter=15,
            batch_size=batch_size,
        )
    elif method == "locmin":
        Preprocess.baseline_correction(
            data_manager=data_manager,
            method="locmin",
            smooth="none",
            span=None,
            upper=False,
            width=11,
            batch_size=batch_size,
        )
    else:
        Preprocess.baseline_correction(
            data_manager=data_manager,
            method="snip",
            m=50,
            decreasing=True,
            baseline_scale=1.0,
            batch_size=batch_size,
        )


def validate_baseline_data(data_manager: MSDataManagerImzML, method="asls", batch_size: int | None = None):
    if method == "asls":
        effective_batch = BATCH_SIZE if batch_size is None else batch_size
        corrected_dm = Preprocess.baseline_correction(
            data_manager=data_manager,
            method="asls",
            lam=1e7,
            p=0.01,
            niter=15,
            batch_size=effective_batch,
        )
    elif method == "locmin":
        effective_batch = BATCH_SIZE if batch_size is None else batch_size
        corrected_dm = Preprocess.baseline_correction(
            data_manager=data_manager,
            method="locmin",
            smooth="none",
            span=None,
            upper=False,
            width=11,
            batch_size=effective_batch,
        )
    else:
        effective_batch = BATCH_SIZE if batch_size is None else batch_size
        corrected_dm = Preprocess.baseline_correction(
            data_manager=data_manager,
            method="snip",
            m=50,
            decreasing=True,
            baseline_scale=1.0,
            batch_size=effective_batch,
        )

    for spectrum, corrected in zip(data_manager.ms, corrected_dm.ms):
        x = spectrum.intensity
        mz_list = spectrum.mz_list
        y = corrected.intensity

        assert x.shape == y.shape
        assert y.shape == mz_list.shape
        assert np.all(y >= 0.0)

class TestAllBaseline:
    @pytest.mark.parametrize("method", [ "locmin","snip"])
    @pytest.mark.benchmark(timer=time.perf_counter)
    def test_all_baseline_methods_replace_data_benchmark(self, benchmark, method, baseline_data_manager):

        logger.info(f"Testing baseline correction with method: {method}")
        benchmark.pedantic(
            partial(run_baseline_correction, method=method, batch_size=BATCH_SIZE),
            args=(baseline_data_manager,),
            rounds=10,
            iterations=1,
            warmup_rounds=0,
        )

        logger.info(f"Finished baseline correction with method: {method}")

    @pytest.mark.parametrize("method", ["locmin", "snip"])
    def test_all_baseline_methods_correctness(self, method, baseline_data_manager):
        validate_baseline_data(baseline_data_manager, method=method, batch_size=BATCH_SIZE)

    def test_asls_correctness_with_subset(self, baseline_data_manager):
        # Build a fresh data manager that only loads a 1/10 spatial window
        subset_ms = MassSpectrumSet()
        dm_subset = MSDataManagerImzML(subset_ms, filepath=baseline_data_manager.filepath)

        dm_subset.extract_metadata()
        height = int(dm_subset.ms.meta.max_count_of_pixels_y)
        width = int(dm_subset.ms.meta.max_count_of_pixels_x)
        subset_height = max(1, height // 10)

        dm_subset.target_locs = [(1, 1), (width, subset_height)]
        dm_subset.load_full_data_from_file()
        for _ in dm_subset.get_batch_generator(batch_size=512):
            pass

        validate_baseline_data(dm_subset, method="asls", batch_size=BATCH_SIZE)