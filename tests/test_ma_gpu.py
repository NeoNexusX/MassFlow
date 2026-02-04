import time
import numpy as np
import pytest
from massflow.module.mass_spectrum_set import MassSpectrumSet
from massflow.module.ms_data_manager_imzml import MSDataManagerImzML
from massflow.tools.logger import get_logger
from massflow.preprocess.dm_pre_fun import Preprocess

# Skip tests if cupy is not available
try:
    import cupy
    CUPY_AVAILABLE = True
except ImportError:
    CUPY_AVAILABLE = False

pytestmark = pytest.mark.filterwarnings("ignore:This process .* is multi-threaded, use of fork():DeprecationWarning")

logger = get_logger("test_noise_reduction_gpu")

@pytest.fixture(scope="module")
def noise_data_manager(data_file_path="Data/other/Example_read/example.imzML") -> MSDataManagerImzML:
    mass_data = MassSpectrumSet()
    dm = MSDataManagerImzML(mass_data, filepath=data_file_path)

    dm.load_full_data_from_file()
    # Pre-load data
    for _ in dm.get_batch_generator(batch_size=512):
        pass
    return dm

def run_dm_noise_reduction_task(
    dm: MSDataManagerImzML,
    method: str = "ma",
    window: int = 11,
    batch_size: int = 256
) -> MSDataManagerImzML:
    denoised_manager = Preprocess.noise_reduction(
        data_manager=dm,
        method=method,
        window=window,
        batch_size=batch_size,
    )
    return denoised_manager

@pytest.mark.skipif(not CUPY_AVAILABLE, reason="CuPy not installed, skipping GPU tests")
class TestNoiseReductionDMGPU:

    @pytest.mark.benchmark(timer=time.perf_counter)
    def test_ma_gpu_speed(
        self,
        benchmark,
        noise_data_manager: MSDataManagerImzML,
    ) -> None:
        """Benchmark GPU implementation of Moving Average."""
        logger.info("Running ma_gpu speed test...")
        denoised_manager = benchmark.pedantic(
            run_dm_noise_reduction_task,
            args=(noise_data_manager, "ma_gpu", 11, 1024),
            rounds=1,
            iterations=10,
            warmup_rounds=1,
        )
        denoised_manager.close()

    # def test_ma_gpu_consistency(self, noise_data_manager: MSDataManagerImzML) -> None:
    #     """Verify that ma_gpu produces results consistent with CPU baseline (ma)."""
    #     dm = noise_data_manager
        
    #     # 1. Run Python baseline (CPU)
    #     logger.info("Running Python baseline (ma) on DM layer...")
    #     dm_cpu = run_dm_noise_reduction_task(
    #         dm, 
    #         method="ma",
    #         window=11, 
    #         batch_size=256
    #     )

    #     # 2. Run GPU target
    #     logger.info("Running GPU target (ma_gpu) on DM layer...")
    #     dm_gpu = run_dm_noise_reduction_task(
    #         dm, 
    #         method="ma_gpu", 
    #         window=11, 
    #         batch_size=256
    #     )

    #     # 3. Compare results
    #     ms_cpu = dm_cpu.ms
    #     ms_gpu = dm_gpu.ms

    #     assert len(ms_cpu) == len(ms_gpu), "Result spectrum count mismatch"

    #     subset_size = min(1000, len(ms_cpu))
    #     logger.info(f"Comparing first {subset_size} spectra between CPU and GPU...")

    #     for i in range(subset_size):
    #         # Compare Intensity
    #         # Note: GPU vs CPU floating point operations might have slight differences
    #         # Relaxing tolerance slightly compared to CPU-CPU tests
    #         np.testing.assert_allclose(
    #             ms_cpu[i].intensity,
    #             ms_gpu[i].intensity,
    #             rtol=1e-5,
    #             atol=1e-5,
    #             err_msg=f"DM GPU consistency failed for method: ma vs ma_gpu at spectrum index {i}"
    #         )
    #         # Compare m/z (should be untouched)
    #         np.testing.assert_array_equal(
    #             ms_cpu[i].mz_list,
    #             ms_gpu[i].mz_list,
    #             err_msg=f"m/z axis mismatch at index {i}"
    #         )

    #     dm_cpu.close()
    #     dm_gpu.close()

    #     logger.info("Consistency check passed for ma vs ma_gpu")