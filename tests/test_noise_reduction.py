import time
import numpy as np
import pytest
from massflow.module.mass_spectrum_set import MassSpectrumSet
from massflow.module.ms_data_manager_imzml import MSDataManagerImzML
from massflow.tools.logger import get_logger
from massflow.preprocess.dm_pre_fun import Preprocess
pytestmark = pytest.mark.filterwarnings("ignore:This process .* is multi-threaded, use of fork():DeprecationWarning")

logger = get_logger("test_noise_reduction")

@pytest.fixture(scope="module")
def noise_data_manager(data_file_path="Data/other/Example_read/example.imzML") -> MSDataManagerImzML:
    """Fixture providing MSDataManagerImzML instance with fully initialized spectra for noise reduction tests."""
    mass_data = MassSpectrumSet()
    dm = MSDataManagerImzML(mass_data, filepath=data_file_path)
    dm.load_full_data_from_file()
    for _ in dm.get_batch_generator(batch_size=512):
        pass
    return dm

def run_noise_reduction_task(
    dm: MSDataManagerImzML,
    method: str = "ma",
    window: int = 11,
    polyorder: int = 3,
    batch_size: int = 256,
) -> MSDataManagerImzML:
    """Execute noise reduction preprocessing task with specified parameters, including batch size."""
    denoised_manager = Preprocess.noise_reduction(
        data_manager=dm,
        method=method,
        window=window,
        polyorder=polyorder,
        batch_size=batch_size,
    )
    return denoised_manager

def validate_noise_reduction_result(
    result_manager: MSDataManagerImzML,
    original_manager: MSDataManagerImzML,
) -> None:
    """Validate basic correctness of noise reduction results: same spectrum count, unchanged mz axis,
    consistent intensity shape, and all finite values."""
    assert len(result_manager.ms) == len(original_manager.ms), "Spectrum count mismatch after noise reduction"

    for i, (src, dst) in enumerate(zip(original_manager.ms, result_manager.ms)):
        assert dst.mz_list is not None and src.mz_list is not None, f"Spectrum {i} has None mz_list"
        assert np.array_equal(dst.mz_list, src.mz_list), f"Spectrum {i} mz_list changed after noise reduction"

        dst_intensity = dst.intensity
        src_intensity = src.intensity
        assert dst_intensity is not None and src_intensity is not None, f"Spectrum {i} has None intensity"
        assert dst_intensity.shape == src_intensity.shape, f"Spectrum {i} intensity shape mismatch after noise reduction"
        assert np.all(np.isfinite(dst_intensity)), f"Spectrum {i} contains NaN or Inf after noise reduction"

class TestNoiseReductionAPI:
    """Test suite for noise reduction API functionality and performance."""

    @pytest.mark.parametrize("method", ["wavelet"])
    @pytest.mark.benchmark(timer=time.perf_counter)
    def test_noise_reduction_benchmark(
        self,
        benchmark,
        noise_data_manager: MSDataManagerImzML,
        method: str,
    ) -> None:
        """Performance test: benchmark Preprocess.noise_reduction execution time."""
        dm_denoised = benchmark.pedantic(
            run_noise_reduction_task,
            args=(noise_data_manager, method, 11, 3, 512),
            rounds=1,
            iterations=1,
            warmup_rounds=0,
        )
        dm_denoised.close()

    @pytest.mark.parametrize("method", ["wavelet"])
    def test_noise_reduction_validation(
        self,
        noise_data_manager: MSDataManagerImzML,
        method: str,
    ) -> None:
        """Functional test: validate correctness of Preprocess.noise_reduction output."""
        dm_denoised = run_noise_reduction_task(
            noise_data_manager,
            method=method,
            window=11,
            polyorder=3,
        )

        validate_noise_reduction_result(dm_denoised, noise_data_manager)
        dm_denoised.close()