import time
import numpy as np
import pytest
from massflow.module.mass_spectrum_set import MassSpectrumSet
from massflow.module.ms_data_manager_imzml import MSDataManagerImzML
from massflow.tools.logger import get_logger
from massflow.preprocess.dm_pre_fun import Preprocess
from massflow.preprocess.spectrum_pre_fun import SpectrumPreprocess
from massflow.preprocess.numba.noise_reduction_numba import (
    smooth_signal_savgol_numba,
    smooth_ns_signal_ma_numba,
    smooth_ns_signal_gaussian_numba,
    smooth_ns_signal_bi_numba,
)
pytestmark = pytest.mark.filterwarnings("ignore:This process .* is multi-threaded, use of fork():DeprecationWarning")

logger = get_logger("test_noise_reduction_numba")


@pytest.fixture(scope="module")
def noise_data_manager(data_file_path="data/example.imzML") -> MSDataManagerImzML:

    mass_data = MassSpectrumSet()
    dm = MSDataManagerImzML(mass_data, filepath=data_file_path)

    dm.extract_metadata()
    height = int(dm.ms.meta.max_count_of_pixels_y)
    width = int(dm.ms.meta.max_count_of_pixels_x)
    subset_height = max(1, height // 10)
    dm.target_locs = [(1, 1), (width, subset_height)]

    dm.load_full_data_from_file()
    for _ in dm.get_batch_generator(batch_size=512):
        pass
    return dm

def run_dm_noise_reduction_task(
    dm: MSDataManagerImzML,
    method: str = "savgol_numba",
    window: int = 11,
    polyorder: int = 3,
    batch_size: int = 256
) -> MSDataManagerImzML:
    denoised_manager = Preprocess.noise_reduction(
        data_manager=dm,
        method=method,
        window=window,
        polyorder=polyorder,
        batch_size=batch_size,
    )
    return denoised_manager

def run_noise_reduction_loop(
    data: MassSpectrumSet,
    method: str = "savgol",
    window: int = 11,
) -> None:
    for i in range(len(data)):
        data[i] = SpectrumPreprocess.noise_reduction_spectrum(
            data=data[i],
            method=method,
            window=window,
        )

class TestNoiseReductionDMNumba:

    @pytest.mark.parametrize(
        "method,backend",
        [
            # ("savgol", "python"),
            # ("savgol", "numba"),
            ("ma_ns", "python"),
            ("ma_ns", "numba"),
            ("gaussian_ns", "python"),
            ("gaussian_ns", "numba"),
            ("bi_ns", "python"),
            ("bi_ns", "numba"),
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
        else:
            dm_method = f"{method}_numba"

        denoised_manager = benchmark.pedantic(
            run_dm_noise_reduction_task,
            args=(noise_data_manager, dm_method, 11, 3, 256),
            rounds=1,
            iterations=1,
            warmup_rounds=0,
        )
        denoised_manager.close()

    @pytest.mark.parametrize(
        "method,numba_kind",
        [
            # ("savgol", "savgol"),
            ("ma_ns", "ma_ns"),
            ("gaussian_ns", "gaussian_ns"),
            ("bi_ns", "bi_ns"),
        ],
    )
    def test_numba_consistency(self, noise_data_manager: MSDataManagerImzML, method: str, numba_kind: str) -> None:
        loaded_mass_data = noise_data_manager.ms
        subset_size = min(1000, len(loaded_mass_data))
        data_for_original = MassSpectrumSet()
        subset_spectra: list = []
        for i in range(subset_size):
            sp = loaded_mass_data[i]
            data_for_original.add_spectrum(sp)
            subset_spectra.append(sp)

        logger.info(
            f"Checking Numba consistency for method={method} (SciPy vs Numba), subset_size={subset_size}"
        )

        # Baseline: Python/SciPy implementation via SpectrumPreprocess
        run_noise_reduction_loop(data_for_original, method=method, window=11)

        if numba_kind == "savgol":
            intensities = np.array(
                [s.intensity for s in subset_spectra],
                dtype=np.float64,
            )
            smoothed_numba = smooth_signal_savgol_numba(
                intensities,
                window=11,
                polyorder=3,
            )

            for i in range(subset_size):
                np.testing.assert_allclose(
                    data_for_original[i].intensity,
                    smoothed_numba[i],
                    rtol=1e-5,
                    atol=1e-5,
                    err_msg=(
                        f"2D Numba consistency failed for method: {method} "
                        f"at spectrum index {i}"
                    ),
                )
        else:
            # NS-series methods: check Numba consistency by calling each method's Numba implementation per spectrum
            for i in range(subset_size):
                sp = subset_spectra[i]
                intensity = sp.intensity
                index = sp.mz_list

                if numba_kind == "ma_ns":
                    numba_intensity = smooth_ns_signal_ma_numba(
                        intensity,
                        index=index,
                        k=11,
                        p=2,
                    )
                elif numba_kind == "gaussian_ns":
                    numba_intensity = smooth_ns_signal_gaussian_numba(
                        intensity,
                        index=index,
                        k=11,
                        p=2,
                        sd=None,
                    )
                else:  # bi_ns
                    numba_intensity = smooth_ns_signal_bi_numba(
                        intensity,
                        index=index,
                        k=11,
                        p=2,
                        sd_dist=None,
                        sd_intensity=None,
                    )

                np.testing.assert_allclose(
                    data_for_original[i].intensity,
                    numba_intensity,
                    rtol=1e-5,
                    atol=1e-5,
                    err_msg=(
                        f"Numba NS consistency failed for method: {method} "
                        f"at spectrum index {i}"
                    ),
                )

        logger.info(
            f"Consistency check passed for method={method} (Numba vs SciPy) "
            f"on first {subset_size} spectra"
        )
