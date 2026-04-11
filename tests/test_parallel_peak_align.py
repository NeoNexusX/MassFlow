import time
from typing import Optional

import pytest
import numpy as np
from numba import set_num_threads

from massflow.tools.logger import get_logger
from massflow.module import MassSpectrumSet
from massflow.data_manager import MSDataManagerImzML
from massflow.preprocess.batch_pre_fun import BatchPreprocess
from massflow.preprocess.helper.peak_align_helper import compute_reference
from massflow.preprocess.helper.peak_align_parallel import compute_reference_parallel, align_spectra_parallel

logger = get_logger("test_peak_alignment")

# FILEPATH = r"/Users/dre/Desktop/data/20170105_ADP_larvae_2-2_NEDC001_400x300_20x20/20170105_ADP_larvae_2-2_NEDC001_400x300_20x20.imzML"
FILEPATH = r"/Users/dre/Desktop/data/ME_Drosophila_Sagittal_Head_DAN_AIF/2023-10-30_ME_Drosophila_Sagittal_Head_DAN_AIF_600-1200_30NCE_100-400_34at_3umss_133x233.imzML"
ROUND = 2

def run_compute_reference_series(dm, units, binfun, binratio, tolerance):
    """Run compute_reference in series mode"""
    return compute_reference(dm, units=units, binfun=binfun, binratio=binratio, tolerance=tolerance)

def run_compute_reference_parallel(dm, units, binfun, binratio, tolerance, numba_max_threads):
    """Run compute_reference in parallel mode"""
    set_num_threads(numba_max_threads)
    return compute_reference_parallel(dm, units=units, binfun=binfun, binratio=binratio, tolerance=tolerance)

def run_alignment_series(
    dm: MSDataManagerImzML,
    units: str,
    tolerance: float,
    reference: np.ndarray,
) -> None:
    """Run peak alignment in series mode"""

    for batch in dm.batch_generator(batch_size=256):
        _ = BatchPreprocess.peak_align_batch(
            batch,
            reference=reference,
            tolerance=tolerance,
            units=units,
        )

def run_alignment_parallel(
    dm: MSDataManagerImzML,
    units: str,
    tolerance: float,
    reference: np.ndarray,
    numba_max_threads: int = 4,
) -> None:
    """Run peak alignment in parallel mode"""
    set_num_threads(numba_max_threads)

    for mz_data, intensity_flat, lengths, _ in dm.flat_generator(batch_size=256):
        _ = align_spectra_parallel(
            mz_data,
            intensity_flat,
            lengths=lengths,
            reference=reference,
            tolerance=tolerance,
            units=units
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
    @pytest.mark.parametrize("binfun", ["median"])
    @pytest.mark.parametrize("binratio", [2])
    @pytest.mark.parametrize("mode, numba_max_threads", [("series", 1), ("parallel", 1), ("parallel", 2), ("parallel", 4)])
    @pytest.mark.benchmark(timer=time.perf_counter)
    def test_compute_reference(
        self,
        benchmark,
        data_manager: MSDataManagerImzML,
        units: str,
        tolerance: Optional[float],
        binfun: str,
        binratio: int,
        mode: str,
        numba_max_threads: int,
    ) -> None:
        """
        Performance Test: Run peak alignment benchmark.
        """
        if mode == "series":
            _, _ = benchmark.pedantic(
                run_compute_reference_series,
                args=(data_manager, units, binfun, binratio, tolerance),
                rounds=ROUND,
                iterations=1,
                warmup_rounds=1,
            )
        else:
            _, _ = benchmark.pedantic(
                run_compute_reference_parallel,
                args=(data_manager, units, binfun, binratio, tolerance, numba_max_threads),
                rounds=ROUND,
                iterations=1,
                warmup_rounds=1,
            )

    @pytest.mark.parametrize("units, tolerance", [("ppm", None)])
    @pytest.mark.parametrize("binfun", ["median"])
    @pytest.mark.parametrize("binratio", [2])
    @pytest.mark.parametrize("mode, numba_max_threads", [("series", 1), ("parallel", 1), ("parallel", 2), ("parallel", 4)])
    @pytest.mark.benchmark(timer=time.perf_counter)
    def test_peak_align(
        self,
        benchmark,
        data_manager: MSDataManagerImzML,
        units: str,
        tolerance: Optional[float],
        binfun: str,
        binratio: int,
        numba_max_threads: int,
        mode: str,
    ) -> None:
        """
        Performance Test: Run parallel peak alignment benchmark.
        """
        if mode == "series":
            reference, tolerance = run_compute_reference_series(data_manager, units, binfun, binratio, tolerance)
            _ = benchmark.pedantic(
                run_alignment_series,
                args=(data_manager, units, tolerance, reference),
                rounds=ROUND,
                iterations=1,
                warmup_rounds=1,
            )
        else:
            reference, tolerance = run_compute_reference_parallel(data_manager, units, binfun, binratio, tolerance, numba_max_threads)
            _ = benchmark.pedantic(
                run_alignment_parallel,
                args=(data_manager, units, tolerance, reference, numba_max_threads),
                rounds=ROUND,
                iterations=1,
                warmup_rounds=1,
            )

    def test_reference_rightness(self, data_manager: MSDataManagerImzML):
        """Test that the reference spectrum is computed correctly"""
        reference, tolerance = compute_reference(data_manager, units="ppm", binfun="median", binratio=2, tolerance=None)
        parallel_reference, parallel_tolerance = compute_reference_parallel(data_manager, units="ppm", binfun="median", binratio=2, tolerance=None)
        assert np.array_equal(reference, parallel_reference), "Reference arrays are not exactly equal!"
        assert tolerance == parallel_tolerance, "Tolerances are not exactly equal!"
