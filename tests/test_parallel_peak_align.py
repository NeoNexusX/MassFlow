import time
from typing import Optional
import pytest
import numpy as np
from massflow.tools.logger import get_logger
from massflow.module import MassSpectrumSet
from massflow.data_manager import MSDataManagerImzML
from massflow.preprocess import Preprocessor
from massflow.preprocess.helper.peak_align_helper import compute_reference_parallel, align_spectra_parallel
from massflow.r_preprocess import set_default_r_home

logger = get_logger("test_peak_alignment")

FILEPATH = r"/Users/dre/Desktop/data/20170105_ADP_larvae_2-2_NEDC001_400x300_20x20.imzML"
MZ_RANGES = [(10, 1000)]
TOLERANCE = 1e-7
ROUND = 2

set_default_r_home("/Library/Frameworks/R.framework/Resources")  # set R_HOME path here

def run_alignment_series(
    dm: MSDataManagerImzML,
    units: str,
    binfun: str,
    binratio: int,
    tolerance: Optional[float] = None,
    backend_method: str = "python",
) -> MSDataManagerImzML:

    align_manager = (
        Preprocessor(dm)
        .peak_align(
            reference=None,
            units=units,
            tolerance=tolerance,
            binfun=binfun,
            binratio=binratio,
            backend=backend_method,
        )
        .start()
    )
    return align_manager

def run_alignment_parallel(
    dm: MSDataManagerImzML,
    units: str,
    binfun: str,
    binratio: int,
    tolerance: Optional[float] = None,
) -> MSDataManagerImzML:

    reference, tolerance = compute_reference_parallel(
        dm, units=units, binfun=binfun, binratio=binratio, tolerance=tolerance
    )

    aligned_ms = MassSpectrumSet()
    aligned_data_manager = MSDataManagerImzML(aligned_ms, temp_dir="./temp_parallel_align")
    aligned_data_manager.copy_meta(dm)
    writer = aligned_data_manager.writer

    for mz_data, intensity_matrix, lengths, coords in dm.matrix_generator(batch_size=256):
        aligned_matrix = align_spectra_parallel(
            mz_data,
            intensity_matrix,
            lengths=lengths,
            reference=reference,
            tolerance=tolerance,
            units=units,
        )
        aligned_data_manager.swap_matrix_data_out2disk(reference, aligned_matrix, lengths, coords)

    aligned_data_manager.close_writer()
    aligned_data_manager.load_head_data()  # Load head data for validation

    return aligned_data_manager

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
    @pytest.fixture(scope="class")
    def data_manager(self) -> MSDataManagerImzML:
        mass_data = MassSpectrumSet()
        dm = MSDataManagerImzML(mass_data, filepath=FILEPATH)
        dm.load_head_data()
        return dm

    @pytest.mark.parametrize("units, tolerance", [("ppm", None)])
    @pytest.mark.parametrize("binfun", ["median", "max"])
    @pytest.mark.parametrize("binratio", [2])
    @pytest.mark.parametrize("backend_method", ["python"])
    @pytest.mark.benchmark(timer=time.perf_counter)
    def test_peak_align_series_benchmark(
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
            run_alignment_series,
            args=(data_manager, units, binfun, binratio, tolerance, backend_method),
            rounds=ROUND,
            iterations=1,
            warmup_rounds=0,
        )

        dm_align.close()

    @pytest.mark.parametrize("units, tolerance", [("ppm", None)])
    @pytest.mark.parametrize("binfun", ["median", "max"])
    @pytest.mark.parametrize("binratio", [2])
    @pytest.mark.benchmark(timer=time.perf_counter)
    def test_peak_align_parallel_benchmark(
        self,
        benchmark,
        data_manager: MSDataManagerImzML,
        units: str,
        tolerance: Optional[float],
        binfun: str,
        binratio: int,
    ) -> None:
        """
        Performance Test: Run parallel peak alignment benchmark.
        """
        dm_align = benchmark.pedantic(
            run_alignment_parallel,
            args=(data_manager, units, binfun, binratio, tolerance),
            rounds=ROUND,
            iterations=1,
            warmup_rounds=0,
        )

        dm_align.close()

    # @pytest.mark.parametrize("units, tolerance", [("ppm", None)])
    # @pytest.mark.parametrize("binfun", ["median"])
    # @pytest.mark.parametrize("binratio", [2])
    # @pytest.mark.parametrize("mz_range", MZ_RANGES)
    # def test_peak_alignment_validation(
    #     self,
    #     data_manager: MSDataManagerImzML,
    #     units: str,
    #     tolerance: Optional[float],
    #     binfun: str,
    #     binratio: int,
    #     mz_range: tuple[float, float],
    # ) -> None:
    #     """
    #     Functional Test: Verify the correctness of peak alignment results.
    #     """

    #     series_align_manager = run_alignment_series(
    #         data_manager, units, binfun, binratio, tolerance
    #     )

    #     parallel_align_manager = run_alignment_parallel(
    #         data_manager, units, binfun, binratio, tolerance
    #     )

    #     compare_align_result(
    #         py_manager=series_align_manager,
    #         r_manager=parallel_align_manager,
    #         mz_range=mz_range
    #     )

    #     series_align_manager.close()
    #     parallel_align_manager.close()
