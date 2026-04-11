import time
import pytest
import numpy as np
from massflow.data_manager import MSDataManagerImzML
from massflow.preprocess.helper.est_noise_helper import estimator
from massflow.tools.logger import get_logger

logger = get_logger("test_est")

ROUNDS = 5
FLAT_EST_DENOISE_METHODS = ["gaussian_numba"]
MAX_SNR_PIXELS = 10


def _estimate_flat_from_flat_batches(
    flat_batches,
    denoise_method: str,
    method: str = "sd",
):
    pixel_counter = 0
    for batch_idx, (intensity_flat, lengths) in enumerate(flat_batches):
        noise_flat = estimator(
            intensity=intensity_flat,
            indexes=None,
            lengths=lengths,
            method=method,
            nbins=1,
            denoise_method=denoise_method,
        )
        offset = 0
        for local_pixel_idx, valid_len in enumerate(lengths):
            end = offset + int(valid_len)
            intensity_slice = intensity_flat[offset:end]
            noise_slice = noise_flat[offset:end]

            signal_level = float(np.percentile(intensity_slice, 95))
            noise_mean = float(np.mean(noise_slice))
            snr = float(signal_level / noise_mean) if noise_mean > 0 else float("inf")

            logger.info(
                "Estimator[%s] batch=%d pixel=%d(local=%d) signal95=%.6f noise=%.6f snr=%.6f",
                denoise_method,
                batch_idx,
                pixel_counter,
                local_pixel_idx,
                signal_level,
                noise_mean,
                snr,
            )

            pixel_counter += 1
            offset = end
            assert snr >= 0, f"SNR should be non-negative, got {snr} for pixel {pixel_counter} in batch {batch_idx}"

            if pixel_counter >= MAX_SNR_PIXELS:
                return


class TestNoiseEstimationAPI:
    """Noise estimation benchmarks focused on flat + numba modes."""

    @pytest.fixture(scope="module", params=["data/example.imzML"])
    def flat_caches(self, request):
        data_file_path = request.param
        dm = MSDataManagerImzML(filepath=data_file_path)
        dm.load_head_data()

        caches = []
        for _, intensity_flat, lengths, _ in dm.flat_generator(
            batch_size=4096,
            include_mz=False,
            max_threads=16,
        ):
            caches.append((intensity_flat, lengths))

        return caches

    @pytest.mark.benchmark(timer=time.perf_counter)
    @pytest.mark.parametrize("denoise_method", FLAT_EST_DENOISE_METHODS)
    def test_est_flat_numba_speed(self, benchmark, denoise_method, flat_caches):
        logger.info("Benchmarking flat estimator denoise_method=%s", denoise_method)

        benchmark.pedantic(
            _estimate_flat_from_flat_batches,
            args=(flat_caches, denoise_method, "sd"),
            rounds=ROUNDS,
            iterations=1,
            warmup_rounds=1,
        )
