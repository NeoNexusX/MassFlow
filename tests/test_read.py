import time
import pytest

from massflow.module.ms_data_manager_imzml import MSDataManagerImzML
from massflow.module.ms_module import MS
from massflow.logger import get_logger
logger = get_logger("test_read")


def data_reading_benchmark(data_file_path="data/example.imzML"):
    mass_data = MS()
    with MSDataManagerImzML(mass_data, filepath=data_file_path) as ms_md:
        ms_md.load_full_data_from_file()
        ms_md.inspect_data()
        for spectrum_data in mass_data:
            _ = spectrum_data.intensity

class TestDataReading:

    @pytest.mark.parametrize("data_file_path", ["data/example.imzML"])
    @pytest.mark.benchmark(timer=time.perf_counter)
    def test_data_reading_benchmark(self, benchmark,data_file_path):
        logger.info("Starting data reading benchmark")
        benchmark.pedantic(data_reading_benchmark,
                           args=(data_file_path,),
                           rounds=10,
                           iterations=1,
                           warmup_rounds=1)
        logger.info("Finished data reading benchmark")