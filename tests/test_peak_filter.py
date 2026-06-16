"""Test peak filtering functionality."""

import os
import time

import pytest

from massflow.data_manager import MSDataManagerImzML
from massflow.preprocess.preprocessor import Preprocessor
from massflow.preprocess.helper.stat_filter_helper import filter_peaks_by_stats
from massflow.tools.logger import get_logger

logger = get_logger("massflow.tests.test_peak_filter")

FILE_PATH = '/Users/dre/Desktop/data/160TopL,130TopR,150BottomL,140BottomR-profile/160TopL,130TopR,150BottomL,140BottomR-profile.imzML'
ROUND = 1

def run_stat_filter(data_manager: MSDataManagerImzML) -> MSDataManagerImzML:
    """Run the statistical filter on the data manager."""

    return filter_peaks_by_stats(data_manager, min_frequency=0.05, intensity_quantile=0.05)


def run_isotope_filter(data_manager: MSDataManagerImzML) -> MSDataManagerImzML:
    """Run the isotope filter on the data manager."""

    return Preprocessor(data_manager).filter_peaks_by_isotopes().start() # type: ignore


class TestPeakFilter:
    """Test the peak filtering functionality."""

    @pytest.fixture(scope="class")
    def data_manager(self) -> MSDataManagerImzML:
        """Fixture to load the data manager."""
        if not os.path.exists(FILE_PATH):
            pytest.skip(f"Missing test data: {FILE_PATH}")
        dm = MSDataManagerImzML(filepath = FILE_PATH)
        dm.load_head_data()
        return dm

    @pytest.mark.benchmark(timer=time.perf_counter)
    def test_stat_filter(self, benchmark, data_manager):
        """Test the statistical filter."""
        if data_manager.ms.meta.continuous is not True:
            continuous_dm = Preprocessor(data_manager).peak_align().start() # type: ignore
        else:
            continuous_dm = data_manager

        filtered_dm = benchmark.pedantic(
            run_stat_filter,
            args=(continuous_dm,),
            rounds=ROUND,
            iterations=1,
            warmup_rounds=0,
        )

        continuous_dm.close()
        filtered_dm.close()

    @pytest.mark.benchmark(timer=time.perf_counter)
    def test_isotope_filter(self, benchmark, data_manager):
        """Test the isotope filter."""
        if data_manager.ms.meta.centroid_spectrum is not True:
            centroid_dm = Preprocessor(data_manager).peak_pick().start() # type: ignore
        else:
            centroid_dm = data_manager

        filtered_dm = benchmark.pedantic(
            run_isotope_filter,
            args=(centroid_dm,),
            rounds=ROUND,
            iterations=1,
            warmup_rounds=0,
        )

        centroid_dm.close()
        filtered_dm.close()
