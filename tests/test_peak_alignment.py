import time
from typing import Optional
import pytest
import numpy as np
from massflow.tools.logger import get_logger
from massflow.module.mass_spectrum_set import MassSpectrumSet
from massflow.module.ms_data_manager_imzml import MSDataManagerImzML
from massflow.preprocess.dm_pre_fun import Preprocess
from massflow.r_preprocess import set_default_r_home

logger = get_logger("test_peak_alignment")

FILEPATH = r"picked.imzML"
MZ_RANGES = [(500, 600)]
TOLERANCE = 0.10
ROUND = 1

set_default_r_home("r_home_path_here")  # set R_HOME path here

@pytest.fixture(scope="module")
def data_manager(filepath=FILEPATH) -> MSDataManagerImzML:
    mass_data = MassSpectrumSet()
    dm = MSDataManagerImzML(mass_data, filepath=filepath)
    dm.load_full_data_from_file()
    return dm

def run_alignment_task(
    dm: MSDataManagerImzML,
    units: str,
    binfun: str,
    binratio: int,
    tolerance: Optional[float] = None,
    backend_method: str = "python",
) -> MSDataManagerImzML:

    align_manager = Preprocess.peak_align(
        data_manager=dm,
        ref=None,
        units=units,
        tolerance=tolerance,
        binfun=binfun,
        binratio=binratio,
        backend=backend_method,
    )
    return align_manager

def validate_align_result(
    result_manager: MSDataManagerImzML, original_manager: MSDataManagerImzML
):
    """Validate the rationality of alignment results"""

    assert len(result_manager.ms) == len(original_manager.ms), "Spectrum count mismatch"

    ref_mz = result_manager.ms[0].mz_list
    assert ref_mz is not None, "Reference m/z axis is None"
    for i, spectrum in enumerate(result_manager.ms):
        assert np.array_equal(
            spectrum.mz_list, ref_mz
        ), f"Spectrum {i} has different m/z axis"

    for i, spectrum in enumerate(result_manager.ms):
        intensity = spectrum.intensity
        assert intensity is not None, f"Spectrum {i} has None intensity"
        assert intensity.shape == ref_mz.shape, f"Spectrum {i} intensity shape mismatch"
        assert np.all(
            np.isfinite(intensity)
        ), f"Spectrum {i} contains NaN or Inf values"

        # Free memory to prevent OOM
        spectrum.clear_data()

    logger.info(
        f"Validation passed: {len(result_manager.ms)} spectra aligned to {len(ref_mz)} peaks"
    )

def compare_align_result(
    py_manager: MSDataManagerImzML,
    r_manager: MSDataManagerImzML,
    mz_range,
    tolerance=TOLERANCE,
):
    """Compare peak alignment results between Python and R implementations to verify correctness"""
    py_ms = py_manager.ms
    r_ms = r_manager.ms

    py_ref = py_ms[0].mz_list
    r_ref = r_ms[0].mz_list

    mz_start, mz_end = mz_range

    py_count = np.count_nonzero((py_ref >= mz_start) & (py_ref <= mz_end))
    r_count = np.count_nonzero((r_ref >= mz_start) & (r_ref <= mz_end))

    py_sum_intensity = 0.0
    py_mask = (py_ref >= mz_start) & (py_ref <= mz_end)
    for spectrum in py_ms:
        intensity = spectrum.intensity
        if intensity is not None and len(intensity) > 0:
            py_sum_intensity += np.sum(intensity[py_mask])
        # Free memory
        spectrum.clear_data()

    r_sum_intensity = 0.0
    r_mask = (r_ref >= mz_start) & (r_ref <= mz_end)
    for spectrum in r_ms:
        intensity = spectrum.intensity
        if intensity is not None and len(intensity) > 0:
            r_sum_intensity += np.sum(intensity[r_mask])
        # Free memory
        spectrum.clear_data()

    # Compare peak counts
    count_rel_diff = abs(py_count - r_count) / max(r_count, 1)
    if count_rel_diff > tolerance:
        pytest.fail(
            f"Peak count in {mz_range}: Python={py_count}, cardinal={r_count}\n"
            f"Relative diff: {count_rel_diff:.2%} (tol={tolerance:.2%})"
        )

    # Compare intensity sums
    intensity_rel_diff = abs(py_sum_intensity - r_sum_intensity) / max(r_sum_intensity, 1e-10)
    if intensity_rel_diff > tolerance:
        pytest.fail(
            f"Total intensity sum in {mz_range}: Python={py_sum_intensity:.2e}, cardinal={r_sum_intensity:.2e}\n"
            f"Relative diff: {intensity_rel_diff:.2%} (tol={tolerance:.2%})"
        )

    logger.info(
        f"Peak count check passed for {mz_range}: "
        f"Python={py_count}, cardinal={r_count}, "
        f"count_diff={count_rel_diff:.2%}"
    )
    logger.info(
        f"Intensity sum check passed for {mz_range}: "
        f"Python={py_sum_intensity:.2e}, cardinal={r_sum_intensity:.2e}, "
        f"intensity_diff={intensity_rel_diff:.2%}"
    )

class TestPeakAlignment:
    """test peak alignment functionality"""

    @pytest.mark.parametrize("units, tolerance", [("ppm", None)])
    @pytest.mark.parametrize("binfun", ["median"])
    @pytest.mark.parametrize("binratio", [2])
    @pytest.mark.parametrize("backend_method", ["python", "cardinal"])
    @pytest.mark.benchmark(timer=time.perf_counter)
    def test_peak_alignment_benchmark(
        self,
        benchmark,
        data_manager: MSDataManagerImzML,
        units: str,
        tolerance: Optional[float],
        binfun: str,
        binratio: int,
        backend_method: str,
    ) -> None:
        """
        Performance Test: Run peak alignment benchmark.
        """
        dm_align = benchmark.pedantic(
            run_alignment_task,
            args=(data_manager, units, binfun, binratio, tolerance, backend_method),
            rounds=ROUND,
            iterations=1,
            warmup_rounds=0,
        )

        dm_align.close()

    @pytest.mark.parametrize("units, tolerance", [("ppm", None)])
    @pytest.mark.parametrize("binfun", ["median"])
    @pytest.mark.parametrize("binratio", [2])
    @pytest.mark.parametrize("mz_range", MZ_RANGES)
    def test_peak_alignment_validation(
        self,
        data_manager: MSDataManagerImzML,
        units: str,
        tolerance: Optional[float],
        binfun: str,
        binratio: int,
        mz_range: tuple[float, float],
    ) -> None:
        """
        Functional Test: Verify the correctness of peak alignment results.
        """

        py_align_manager = run_alignment_task(
            data_manager, units, binfun, binratio, tolerance, backend_method="python"
        )

        r_align_manager = run_alignment_task(
            data_manager, units, binfun, binratio, tolerance, backend_method="cardinal"
        )

        validate_align_result(py_align_manager, data_manager)
        validate_align_result(r_align_manager, data_manager)

        compare_align_result(
            py_manager=py_align_manager,
            r_manager=r_align_manager,
            mz_range=mz_range
        )

        py_align_manager.close()
        r_align_manager.close()
