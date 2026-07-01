"""Test peak filtering functionality."""

import os
import time
from pathlib import Path

import pytest

from massflow.data_manager import MSDataManagerImzML
from massflow.preprocess.preprocessor import Preprocessor
from massflow.preprocess.helper.stat_filter_helper import filter_peaks_by_stats
from massflow.tools.logger import get_logger

logger = get_logger("massflow.tests.test_peak_filter")

FILE_PATH = 'data/example.imzML'
ROUND = 1
TEMP_OUTPUT_BASES = (
    "temp/stat_filter",
    "temp/isotope_filter",
    "temp/peak_align_for_stat_filter",
    "temp/peak_pick_for_isotope_filter",
)
TEMP_OUTPUT_SUFFIXES = (".imzML", ".ibd", ".fdata", ".log", ".metadata", ".IBD")


def run_stat_filter(data_manager: MSDataManagerImzML) -> MSDataManagerImzML:
    """Run the statistical filter on the data manager."""

    return filter_peaks_by_stats(data_manager, output_path="temp/stat_filter.imzML", min_frequency=0.05, intensity_quantile=0.05)


def run_isotope_filter(data_manager: MSDataManagerImzML) -> MSDataManagerImzML:
    """Run the isotope filter on the data manager."""

    return Preprocessor(data_manager, output_path="temp/isotope_filter.imzML").filter_peaks_by_isotopes().start() # type: ignore


def del_temp_file():
    """Delete generated test output files if they exist."""
    for base in TEMP_OUTPUT_BASES:
        base_path = Path(base)
        for suffix in TEMP_OUTPUT_SUFFIXES:
            base_path.with_suffix(suffix).unlink(missing_ok=True)

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
        del_temp_file()
        if data_manager.ms.meta.continuous is not True:
            continuous_dm = Preprocessor(data_manager, output_path="temp/peak_align_for_stat_filter.imzML").peak_align().start() # type: ignore
        else:
            continuous_dm = data_manager

        filtered_dm = None
        try:
            filtered_dm = benchmark.pedantic(
                run_stat_filter,
                args=(continuous_dm,),
                rounds=ROUND,
                iterations=1,
                warmup_rounds=0,
            )
        finally:
            if continuous_dm is not data_manager:
                continuous_dm.close()
            if filtered_dm is not None:
                filtered_dm.close()
            del_temp_file()

    @pytest.mark.benchmark(timer=time.perf_counter)
    def test_isotope_filter(self, benchmark, data_manager):
        """Test the isotope filter."""
        del_temp_file()
        if data_manager.ms.meta.centroid_spectrum is not True:
            centroid_dm = Preprocessor(data_manager, output_path="temp/peak_pick_for_isotope_filter.imzML").peak_pick().start() # type: ignore
        else:
            centroid_dm = data_manager

        filtered_dm = None
        try:
            filtered_dm = benchmark.pedantic(
                run_isotope_filter,
                args=(centroid_dm,),
                rounds=ROUND,
                iterations=1,
                warmup_rounds=0,
            )
        finally:
            if centroid_dm is not data_manager:
                centroid_dm.close()
            if filtered_dm is not None:
                filtered_dm.close()
            del_temp_file()
