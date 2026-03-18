import time
from typing import Optional
import numpy as np
import pytest
from massflow.module.mass_spectrum_set import MassSpectrumSet
from massflow.data_manager.ms_data_manager_imzml import MSDataManagerImzML
from massflow.preprocess.dm_pre_fun import Preprocess
from massflow.tools.logger import get_logger

pytestmark = pytest.mark.filterwarnings(
    "ignore:This process .* is multi-threaded, use of fork():DeprecationWarning"
)

logger = get_logger("test_baseline_correction_numba")


@pytest.fixture(scope="module")
def baseline_numba_data_manager(
    data_file_path: str = "Data/other/Example_read/example.imzML",
) -> MSDataManagerImzML:
    mass_data = MassSpectrumSet()
    dm = MSDataManagerImzML(mass_data, filepath=data_file_path)
    dm.load_head_data()

    for _ in dm.batch_generator(batch_size=512):
        pass
    
    logger.info("data pre-load finished!")
    return dm


def run_dm_baseline_task(
    dm: MSDataManagerImzML,
    method: str = "snip",
    batch_size: int = 512,
    numba_max_threads: Optional[int] = None,
) -> MSDataManagerImzML:
    corrected_manager = Preprocess.baseline_correction(
        data_manager=dm,
        method=method,
        m=50,
        decreasing=True,
        baseline_scale=1.0,
        batch_size=batch_size,
        numba_max_threads=numba_max_threads,
    )
    return corrected_manager


class TestBaselineCorrectionDMNumba:
    @pytest.mark.parametrize(
        "method,backend",
        [
            #("snip", "python"),
            ("snip", "numba"),
        ],
    )
    @pytest.mark.benchmark(timer=time.perf_counter)
    def test_snip_numba_speed(
        self,
        benchmark,
        baseline_numba_data_manager: MSDataManagerImzML,
        method: str,
        backend: str,
    ) -> None:
        if backend == "python":
            dm_method = method                # "snip"
            threads = None
        else:
            dm_method = f"{method}_numba"     # "snip_numba"
            threads = 8

        logger.info(
            f"Benchmarking baseline correction with method={dm_method}, "
            f"backend={backend}, threads={threads}"
        )

        corrected_manager = benchmark.pedantic(
            run_dm_baseline_task,
            args=(baseline_numba_data_manager, dm_method, 512, threads),
            rounds=1,
            iterations=3,
            warmup_rounds=0,
        )
        corrected_manager.close()

    def test_snip_numba_consistency(
        self,
        baseline_numba_data_manager: MSDataManagerImzML,
    ) -> None:
        dm = baseline_numba_data_manager

        logger.info("Running baseline SNIP (python) on DM layer for consistency check...")
        dm_python = run_dm_baseline_task(dm, method="snip", batch_size=512)

        logger.info("Running baseline SNIP (numba) on DM layer for consistency check...")
        dm_numba = run_dm_baseline_task(dm, method="snip_numba", batch_size=512)

        ms_python = dm_python.ms
        ms_numba = dm_numba.ms

        assert len(ms_python) == len(ms_numba)

        subset_size = min(1000, len(ms_python))
        logger.info(
            f"Comparing first {subset_size} spectra between python SNIP "
            f"and numba SNIP on DM layer..."
        )

        for i in range(subset_size):
            intensity_python = ms_python[i].intensity
            intensity_numba = ms_numba[i].intensity

            assert intensity_python is not None, f"ms_python[{i}].intensity is None"
            assert intensity_numba is not None, f"ms_numba[{i}].intensity is None"

            np.testing.assert_allclose(
                intensity_python,
                intensity_numba,
                rtol=1e-5,
                atol=1e-5,
                err_msg=(
                    f"SNIP Numba DM consistency failed at spectrum index {i}"
                ),
            )

        dm_python.close()
        dm_numba.close()

        logger.info(
            f"SNIP Numba DM consistency check passed on first {subset_size} spectra."
        )