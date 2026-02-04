import time
from typing import Sequence
import pytest
import numpy as np
from massflow.module.mass_spectrum_set import MassSpectrumSet
from massflow.module.ms_data_manager_imzml import MSDataManagerImzML
from massflow.preprocess.dm_pre_fun import Preprocess
from massflow.r_preprocess import set_default_r_home
from massflow.tools.logger import get_logger

logger = get_logger("test_peak_pick")

FILEPATH = r"./data/example.imzML"
MZ_RANGES = [(500, 600)]
TOLERANCE = 0.10
ROUND = 1

set_default_r_home("r_home_path_here")  # set R_HOME path here

def run_pick_test(
    data: MSDataManagerImzML,
    width: int | Sequence[int] = 2,
    method: str = 'scipy',
    relheight: float = 0.1,
    snr: float = 3.0,
    return_type: str = 'height',
    batch_size: int = 256,
    use_numba: bool = True,
    backend: str = "python"
    ) -> MSDataManagerImzML:
    """Run peak picking test using old or new method."""

    return Preprocess.peak_pick(data_manager=data,
                                width=width,
                                method=method,
                                relheight=relheight,
                                snr=snr,
                                return_type=return_type,
                                batch_size=batch_size,
                                use_numba=use_numba,
                                backend=backend,
                                )

def compare_pick_result(
    py_manager: MSDataManagerImzML,
    r_manager: MSDataManagerImzML,
    mz_range,
    tolerance = TOLERANCE
):
    """Compare peak picking results between python and R implementations."""

    py_ms = py_manager.ms
    r_ms = r_manager.ms

    py_sum_intensity = 0.0
    py_sum_count = 0
    for spectrum in py_ms:
        mz_mask = (spectrum.mz_list >= mz_range[0]) & (spectrum.mz_list <= mz_range[1])
        target_intensity = spectrum.intensity[mz_mask]
        valid_mask = target_intensity > 0
        py_sum_intensity += np.sum(target_intensity[valid_mask])
        py_sum_count += np.count_nonzero(valid_mask)
        spectrum.clear_data()

    r_sum_intensity = 0.0
    r_sum_count = 0
    for spectrum in r_ms:
        mz_mask = (spectrum.mz_list >= mz_range[0]) & (spectrum.mz_list <= mz_range[1])
        target_intensity = spectrum.intensity[mz_mask]
        valid_mask = target_intensity > 0
        r_sum_intensity += np.sum(target_intensity[valid_mask])
        r_sum_count += np.count_nonzero(valid_mask)
        spectrum.clear_data()

    count_rel_diff = abs(py_sum_count - r_sum_count) / max(r_sum_count, 1)
    if count_rel_diff > tolerance:
        pytest.fail(
            f"Peak count in {mz_range}: Python={py_sum_count}, cardinal={r_sum_count}\n"
            f"Relative diff: {count_rel_diff:.2%} (tol={tolerance:.2%})"
        )

    intensity_rel_diff = abs(py_sum_intensity - r_sum_intensity) / max(r_sum_intensity, 1e-10)
    if intensity_rel_diff > tolerance:
        pytest.fail(
            f"Total intensity sum in {mz_range}: Python={py_sum_intensity:.2e}, cardinal={r_sum_intensity:.2e}\n"
        )

    logger.info(
        f"Peak count check passed for {mz_range}: "
        f"Python={py_sum_count}, cardinal={r_sum_count},"
        f"count_diff={count_rel_diff:.2%}"
    )

    logger.info(
        f"Intensity sum check passed for {mz_range}: "
        f"Python={py_sum_intensity:.2e}, cardinal={r_sum_intensity:.2e}, "
        f"intensity_diff={intensity_rel_diff:.2%}"
    )


class TestPeakPick:
    """Test peak picking functions."""

    @pytest.fixture(scope="class")
    def data_manager(self) -> MSDataManagerImzML:
        mass_data = MassSpectrumSet()
        dm = MSDataManagerImzML(mass_data, filepath=FILEPATH)
        dm.load_full_data_from_file()
        return dm

    @pytest.mark.parametrize("width", [2])
    @pytest.mark.parametrize("relheight", [0.012])
    @pytest.mark.parametrize("snr", [2.0])
    @pytest.mark.parametrize("return_type", ["area"])
    @pytest.mark.parametrize("batch_size", [256])
    @pytest.mark.parametrize("use_numba", [True])
    @pytest.mark.parametrize("method, backend", [('origin', 'python'), ('diff', 'cardinal')])
    @pytest.mark.benchmark(timer=time.perf_counter)
    def test_peak_pick_benchmark(
        self,
        benchmark,
        data_manager: MSDataManagerImzML, # pylint: disable=redefined-outer-name
        width: int | Sequence[int],
        method: str,
        relheight: float,
        snr: float,
        return_type: str,
        batch_size: int,
        use_numba: bool,
        backend: str,
    ) -> None:
        """Performance Test: Run peak picking benchmark."""

        dm_picked = benchmark.pedantic(
            run_pick_test,
            args=(
                data_manager,
                width,
                method,
                relheight,
                snr,
                return_type,
                batch_size,
                use_numba,
                backend),
            rounds=ROUND,
            iterations=1,
            warmup_rounds=0)

        dm_picked.close()

    @pytest.mark.parametrize("width", [1])
    @pytest.mark.parametrize("relheight", [0.012])
    @pytest.mark.parametrize("snr", [2.0])
    @pytest.mark.parametrize("return_type", ["height"])
    @pytest.mark.parametrize("py_method", ['origin'])
    @pytest.mark.parametrize("r_method", ['diff'])
    @pytest.mark.parametrize("mz_range", MZ_RANGES)
    def test_peak_pick_validation(
        self,
        data_manager: MSDataManagerImzML, # pylint: disable=redefined-outer-name
        width: int | Sequence[int],
        py_method: str,
        r_method: str,
        relheight: float,
        snr: float,
        return_type: str,
        mz_range: tuple[float, float],
    ) -> None:
        """Functional Test: Verify the correctness of peak picking results."""

        py_picked = run_pick_test(
            data=data_manager,
            width=width,
            method=py_method,
            relheight=relheight,
            snr=snr,
            return_type=return_type,
            backend="python")

        r_picked = run_pick_test(
            data=data_manager,
            width=width,
            method=r_method,
            relheight=relheight,
            snr=snr,
            return_type=return_type,
            backend="cardinal")

        compare_pick_result(
            py_manager=py_picked,
            r_manager=r_picked,
            mz_range=mz_range
        )

        py_picked.close()
        r_picked.close()
