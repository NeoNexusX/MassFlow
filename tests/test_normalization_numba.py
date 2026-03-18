import time
from typing import Optional
import numpy as np
import pytest
from massflow.module.mass_spectrum_set import MassSpectrumSet
from massflow.data_manager.ms_data_manager_imzml import MSDataManagerImzML
from massflow.tools.logger import get_logger
from massflow.preprocess.dm_pre_fun import Preprocess

pytestmark = pytest.mark.filterwarnings("ignore:This process .* is multi-threaded, use of fork():DeprecationWarning")

logger = get_logger("test_normalization_numba")


@pytest.fixture(scope="module")
def norm_data_manager(data_file_path="Data/other/Example_read/example.imzML") -> MSDataManagerImzML:
    """Load example data once per module session."""
    mass_data = MassSpectrumSet()
    dm = MSDataManagerImzML(mass_data, filepath=data_file_path)

    # Note: Using full data load similar to noise reduction test for real-world simulation
    dm.load_head_data()
    # Force generator iteration to ensure data is likely cached/ready if applicable
    for _ in dm.batch_generator(batch_size=512):
        pass
    return dm

def run_dm_normalization_task(
    dm: MSDataManagerImzML,
    method: str = "tic",
    scale_method: str = "none",
    scale: float = 1.0,
    batch_size: int = 256,
    numba_max_threads: Optional[int] = None
) -> MSDataManagerImzML:
    """Helper to run normalization via the high-level DM Preprocess API."""
    normalized_manager = Preprocess.normalization(
        data_manager=dm,
        method=method,
        scale_method=scale_method,
        scale=scale,
        batch_size=batch_size,
        numba_max_threads=numba_max_threads
    )
    return normalized_manager


class TestNormalizationDMNumba:

    @pytest.mark.parametrize(
        "method,backend",
        [
            ("tic", "python"),
            ("tic", "numba"),
            ("rms", "python"),
            ("rms", "numba"),
            ("median", "python"),
            ("median", "numba"),
        ],
    )
    @pytest.mark.benchmark(timer=time.perf_counter)
    def test_norm_numba_speed(
        self,
        benchmark,
        norm_data_manager: MSDataManagerImzML,
        method: str,
        backend: str,
    ) -> None:
        """Benchmark speed difference between Python and Numba backends."""
        if backend == "python":
            dm_method = method
            threads = None
        else:
            # Construct numba method name: e.g. "tic" -> "tic_numba"
            dm_method = f"{method}_numba"
            threads = 4

        # Run benchmark
        normalized_manager = benchmark.pedantic(
            run_dm_normalization_task,
            args=(norm_data_manager, dm_method, "none", 1.0, 512, threads),
            rounds=1, # Adjust based on data size
            iterations=1,
            warmup_rounds=0,
        )
        normalized_manager.close()

    @pytest.mark.parametrize(
        "method",
        ["tic", "rms", "median"],
    )
    def test_numba_consistency(self, norm_data_manager: MSDataManagerImzML, method: str) -> None:
        """Verify that Numba results match Python baseline results."""
        dm = norm_data_manager
        
        # 1. Run Python baseline
        logger.info(f"Running Python baseline for method={method} on DM layer...")
        dm_python = run_dm_normalization_task(
            dm, 
            method=method, 
            batch_size=256
        )

        # 2. Run Numba target
        dm_method_numba = f"{method}_numba"
        logger.info(f"Running Numba target for method={dm_method_numba} on DM layer...")
        dm_numba = run_dm_normalization_task(
            dm, 
            method=dm_method_numba, 
            batch_size=256,
            numba_max_threads=4 
        )

        # 3. Compare results
        ms_python = dm_python.ms
        ms_numba = dm_numba.ms

        assert len(ms_python) == len(ms_numba), "Result spectrum count mismatch"

        subset_size = min(1000, len(ms_python))
        logger.info(f"Comparing first {subset_size} spectra between Python and Numba...")

        for i in range(subset_size):
            # Compare Intensity
            # Note: Floating point arithmetic order might differ slightly between compiled C and Python Sum
            # Using rtol=1e-5 should be safe enough for float64
            intensity_py = ms_python[i].intensity
            intensity_numba = ms_numba[i].intensity
            assert intensity_py is not None, f"Python intensity is None at index {i}"
            assert intensity_numba is not None, f"Numba intensity is None at index {i}"
            np.testing.assert_allclose(
                intensity_py,
                intensity_numba,
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