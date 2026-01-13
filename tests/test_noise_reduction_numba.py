import time
import copy
import numpy as np
import pytest
from massflow.module.mass_spectrum_set import MassSpectrumSet
from massflow.module.ms_data_manager_imzml import MSDataManagerImzML
from massflow.tools.logger import get_logger
from massflow.preprocess.dm_pre_fun import Preprocess
from massflow.preprocess.spectrum_pre_fun import SpectrumPreprocess
from massflow.preprocess.numba.noise_reduction_numba import smooth_signal_savgol_numba
pytestmark = pytest.mark.filterwarnings("ignore:This process .* is multi-threaded, use of fork():DeprecationWarning")

logger = get_logger("test_noise_reduction_numba")


@pytest.fixture(scope="module")
def noise_data_manager(data_file_path="data/example.imzML") -> MSDataManagerImzML:

    mass_data = MassSpectrumSet()
    dm = MSDataManagerImzML(mass_data, filepath=data_file_path)
    dm.load_full_data_from_file()
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

    @pytest.mark.benchmark(timer=time.perf_counter)
    def test_dm_savgol_numba_benchmark(
        self,
        benchmark,
        noise_data_manager: MSDataManagerImzML,
    ) -> None:
        dm_denoised = benchmark.pedantic(
            run_dm_noise_reduction_task,
            args=(noise_data_manager, "savgol_numba", 11, 3, 256),
            rounds=1,
            iterations=1,
            warmup_rounds=0,
        )
        dm_denoised.close()

    @pytest.mark.parametrize("method", ["savgol"])
    def test_numba_consistency(self, noise_data_manager: MSDataManagerImzML, method: str) -> None:
        loaded_mass_data = noise_data_manager.ms
        subset_size = min(20, len(loaded_mass_data))
        data_for_original = MassSpectrumSet()
        subset_spectra: list = []
        for i in range(subset_size):
            sp = copy.deepcopy(loaded_mass_data[i])
            data_for_original.add_spectrum(copy.deepcopy(loaded_mass_data[i]))
            subset_spectra.append(sp)
        numba_method = f"{method}_numba"
        logger.info(
            f"Checking 2D Numba consistency for method={method} vs {numba_method}, subset_size={subset_size}"
        )
        run_noise_reduction_loop(data_for_original, method=method)
        intensities = np.array([s.intensity for s in subset_spectra], dtype=np.float64)
        smoothed_2d = smooth_signal_savgol_numba(intensities, window=11, polyorder=3)
        for i in range(subset_size):
            np.testing.assert_allclose(
                data_for_original[i].intensity,
                smoothed_2d[i],
                rtol=1e-5,
                atol=1e-5,
                err_msg=(
                    f"2D Numba consistency failed for method: {method} "
                    f"at spectrum index {i}"
                ),
            )
        logger.info(
            f"Consistency check passed for method={method} (2D Numba vs SciPy) "
            f"on first {subset_size} spectra"
        )

