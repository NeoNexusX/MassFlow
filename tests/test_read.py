import time
import pytest
import gc
from massflow.module.ms_data_manager_imzml import MSDataManagerImzML
from massflow.module.ms_module import MS
from massflow.logger import get_logger
logger = get_logger("test_read")
ROUND = 10
def data_reading_benchmark(data_file_path="data/example.imzML"):
    mass_data = MS()
    with MSDataManagerImzML(mass_data, filepath=data_file_path) as ms_md:
        ms_md.load_full_data_from_file()
        ms_md.inspect_data()
        for index, spectrum_data in enumerate(mass_data):
            _ = spectrum_data.intensity

def data_batch_reading_benchmark(data_file_path="data/example.imzML", batch_size=1024, max_workers=2):
    mass_data = MS()
    with MSDataManagerImzML(mass_data, filepath=data_file_path) as ms_md:
        ms_md.load_full_data_from_file()
        ms_md.inspect_data()
        batch_generator = ms_md.get_batch_generator(batch_size=batch_size, max_workers=max_workers)
        for batch in batch_generator:
            pass


class TestDataReading:

    @pytest.mark.parametrize("data_file_path", ["data/example.imzML"])
    @pytest.mark.benchmark(timer=time.perf_counter)
    def test_data_reading_benchmark(self, benchmark,data_file_path):
        logger.info("Starting data reading benchmark")
        benchmark.pedantic(data_reading_benchmark,
                           args=(data_file_path,),
                           rounds=ROUND,
                           iterations=1,
                           warmup_rounds=1)
        logger.info("Finished data reading benchmark")

    @pytest.mark.parametrize(
        "data_file_path,batch_size,max_workers",
        [
            ("data/example.imzML", 32, 1),
            ("data/example.imzML", 32, 2),
            ("data/example.imzML", 32, 4),
            ("data/example.imzML", 32, 8),
            ("data/example.imzML", 32, 16),
            ("data/example.imzML", 32, 32),
            ("data/example.imzML", 64, 1),
            ("data/example.imzML", 64, 32),
            ("data/example.imzML", 64, 64),
            ("data/example.imzML", 128, 1),
            ("data/example.imzML", 128, 2),
            ("data/example.imzML", 128, 8),
            ("data/example.imzML", 128, 16),
            ("data/example.imzML", 128, 32),
            ("data/example.imzML", 256, 1),
            ("data/example.imzML", 256, 2),
            ("data/example.imzML", 256, 8),
            ("data/example.imzML", 256, 16),
            ("data/example.imzML", 256, 32),
            ("data/example.imzML", 512, 2),
            ("data/example.imzML", 512, 8),
            ("data/example.imzML", 512, 16),
            ("data/example.imzML", 512, 32),
            ("data/example.imzML", 1024, 2),
            ("data/example.imzML", 2048, 2),
            ("data/example.imzML", 2048, 6),
            ("data/example.imzML", 2048, 16),
            ("data/example.imzML", 2048, 32),
            ("data/example.imzML", 2048, 64),
            ("data/example.imzML", 2048, 64),
            ("data/example.imzML", 4068, 64),
            ("data/example.imzML", 4068, 128),
        ]
    )
    @pytest.mark.benchmark(timer=time.perf_counter)
    def test_batch_reading_benchmark(self, benchmark,data_file_path,batch_size,max_workers):
        logger.info("Starting data batch reading benchmark")
        benchmark.pedantic(data_batch_reading_benchmark,
                           args=(data_file_path,batch_size,max_workers),
                           rounds=ROUND,
                           iterations=1,
                           warmup_rounds=1)
        logger.info("Finished data batch reading benchmark")