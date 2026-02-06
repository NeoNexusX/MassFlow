import time
import numpy as np
import pytest
from massflow.module.mass_spectrum_set import MassSpectrumSet
from massflow.module.ms_data_manager_imzml import MSDataManagerImzML
from massflow.tools.logger import get_logger
from massflow.preprocess.dm_pre_fun import Preprocess

pytestmark = pytest.mark.filterwarnings("ignore:This process .* is multi-threaded, use of fork():DeprecationWarning")

logger = get_logger("test_noise_reduction_numba")


@pytest.fixture(scope="module")
def noise_data_manager(data_file_path="Data/other/Example_read/example.imzML") -> MSDataManagerImzML:

    mass_data = MassSpectrumSet()
    dm = MSDataManagerImzML(mass_data, filepath=data_file_path)

    #dm.extract_metadata()
    # height = int(dm.ms.meta.max_count_of_pixels_y)
    # width = int(dm.ms.meta.max_count_of_pixels_x)
    # subset_height = max(1, height // 10)
    # dm.target_locs = [(1, 1), (width, subset_height)]

    dm.load_full_data_from_file()
    for _ in dm.get_batch_generator(batch_size=512):
        pass
    return dm

def run_dm_noise_reduction_task(
    dm: MSDataManagerImzML,
    method: str = "savgol_numba",
    window: int = 11,
    polyorder: int = 3,
    batch_size: int = 256,
    numba_max_threads: int = None
) -> MSDataManagerImzML:
    denoised_manager = Preprocess.noise_reduction(
        data_manager=dm,
        method=method,
        window=window,
        polyorder=polyorder,
        batch_size=batch_size,
        numba_max_threads=numba_max_threads
    )
    return denoised_manager


class TestNoiseReductionDMNumba:

    @pytest.mark.parametrize(
        "method,backend",
        [
            #("savgol", "python"),
            ("savgol", "numba"),
            #("ma", "python"),
            ("ma", "numba"),
            ("ma_loop", "numba"), # New O(N) loop method
            #("gaussian", "python"),
            ("gaussian", "numba"),
            # ("ma_ns", "python"),
            #("ma_ns", "numba"),
            #("gaussian_ns", "python"),
            #("gaussian_ns", "numba"),
            #("bi_ns", "python"),
            #("bi_ns", "numba"),
        
        ],
    )
    @pytest.mark.benchmark(timer=time.perf_counter)
    def test_ns_numba_speed(
        self,
        benchmark,
        noise_data_manager: MSDataManagerImzML,
        method: str,
        backend: str,
    ) -> None:
        if backend == "python":
            dm_method = method
            threads = None
        else:
            # Map simplified method names to actual backend method names
            if method == "ma":
                dm_method = "ma_numba"
            elif method == "ma_loop":
                dm_method = "ma_loop"
            elif method == "gaussian":
                dm_method = "gaussian_numba"
            else:
                # Default fallback for other methods (like savgol -> savgol_numba)
                dm_method = f"{method}_numba"
            threads = 1

        denoised_manager = benchmark.pedantic(
            run_dm_noise_reduction_task,
            args=(noise_data_manager, dm_method, 11, 3, 512, threads),
            rounds=1,
            iterations=10,
            warmup_rounds=1,
        )
        denoised_manager.close()

    @pytest.mark.parametrize(
        "method",
        [
            "ma_ns",
            "gaussian_ns",
            "bi_ns",
            "ma",
            "gaussian",
            "savgol",
            "ma_loop", 
        ],
    )
    def test_numba_consistency(self, noise_data_manager: MSDataManagerImzML, method: str) -> None:
        dm = noise_data_manager
        
        # 1. Run Python baseline (DM layer)
        # For new methods, map back to original python method name for baseline
        if method == "ma_loop":
            baseline_method = "ma"
        else:
            baseline_method = method

        logger.info(f"Running Python baseline for method={baseline_method} on DM layer...")
        dm_python = run_dm_noise_reduction_task(
            dm, 
            method=baseline_method, # e.g. "ma", "gaussian"
            window=11, 
            polyorder=3, 
            batch_size=256
        )

        # 2. Run Numba target (DM layer)
        # Construct the numba method name
        if method == "ma_loop":
            dm_method_numba = "ma_loop"
        else:
            # e.g. "ma" -> "ma_numba", "ma_ns" -> "ma_ns_numba"
            dm_method_numba = f"{method}_numba"

        logger.info(f"Running Numba target for method={dm_method_numba} on DM layer...")
        dm_numba = run_dm_noise_reduction_task(
            dm, 
            method=dm_method_numba, 
            window=11, 
            polyorder=3, 
            batch_size=256,
            numba_max_threads=4 # Use same threads setting as benchmark
        )

        # 3. Compare results
        ms_python = dm_python.ms
        ms_numba = dm_numba.ms

        assert len(ms_python) == len(ms_numba), "Result spectrum count mismatch"

        subset_size = min(1000, len(ms_python))
        logger.info(f"Comparing first {subset_size} spectra between Python and Numba...")

        for i in range(subset_size):
            # Compare Intensity
            np.testing.assert_allclose(
                ms_python[i].intensity,
                ms_numba[i].intensity,
                rtol=1e-5,
                atol=1e-5,
                err_msg=(
                    f"DM Numba consistency failed for method: {method} vs {dm_method_numba} "
                    f"at spectrum index {i}"
                ),
            )
            # Compare m/z (should be untouched)
            np.testing.assert_array_equal(
                ms_python[i].mz_list,
                ms_numba[i].mz_list,
                err_msg=f"m/z axis mismatch at index {i}"
            )

        dm_python.close()
        dm_numba.close()

        logger.info(
            f"Consistency check passed for method={method} vs {dm_method_numba} "
            f"on first {subset_size} spectra"
        )
