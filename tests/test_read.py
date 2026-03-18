import time
import subprocess
import pytest
from massflow.data_manager.ms_data_manager_imzml import MSDataManagerImzML
from massflow.tools.logger import get_logger

logger = get_logger("test_read")

ROUND = 5

def drop_caches():
    subprocess.run(['sudo', 'purge'], check=True)

def data_reading_benchmark(data_file_path="data/example.imzML"):
    with MSDataManagerImzML(filepath=data_file_path) as ms_md:
        ms_md.load_head_data()
        for spectrum_data in ms_md.ms:
            assert spectrum_data.intensity is not None


def data_batch_reading_full_benchmark(
    data_file_path="data/example.imzML", batch_size=1024, max_workers=2
):
    with MSDataManagerImzML(filepath=data_file_path) as ms_md:
        ms_md.load_head_data()
        batch_generator = ms_md.batch_generator(
            batch_size=batch_size,
            max_threads=max_workers
        )
        for batch in batch_generator:
            assert len(batch[1].mz_list) > 0 #type: ignore
            # ms_md.clear_batch_data_memory(batch)

class TestDataReading:

    """
    Test suite for benchmarking data reading performance of MSDataManagerImzML.
    """

    @pytest.mark.parametrize("data_file_path", ["./data/example.imzML"])
    @pytest.mark.benchmark(timer=time.perf_counter)
    def test_data_reading_benchmark(self, benchmark, data_file_path):
        logger.info("Starting data reading benchmark")
        benchmark.pedantic(
            data_reading_benchmark,
            args=(data_file_path,),
            rounds=ROUND,
            iterations=1,
            warmup_rounds=1,
        )
        logger.info("Finished data reading benchmark")

    @pytest.mark.parametrize(
        "data_file_path,batch_size,max_workers",
        [
            # batch_size = 32
            ("data/example.imzML", 32, 1),
            ("data/example.imzML", 32, 2),
            ("data/example.imzML", 32, 8),
            ("data/example.imzML", 32, 32),
            # batch_size = 64
            ("data/example.imzML", 64, 1),
            ("data/example.imzML", 64, 2),
            ("data/example.imzML", 64, 8),
            ("data/example.imzML", 64, 32),
            ("data/example.imzML", 64, 64),
            # batch_size = 128
            ("data/example.imzML", 128, 1),
            ("data/example.imzML", 128, 2),
            ("data/example.imzML", 128, 8),
            ("data/example.imzML", 128, 32),
            ("data/example.imzML", 128, 64),
            # batch_size = 256
            ("data/example.imzML", 256, 1),
            ("data/example.imzML", 256, 2),
            ("data/example.imzML", 256, 8),
            ("data/example.imzML", 256, 32),
            ("data/example.imzML", 256, 64),
            # # batch_size = 512
            # ("data/example.imzML", 512, 1),
            # ("data/example.imzML", 512, 2),
            # ("data/example.imzML", 512, 8),
            # ("data/example.imzML", 512, 32),
            # ("data/example.imzML", 512, 64),
            # # batch_size = 1024
            # ("data/example.imzML", 1024, 1),
            # ("data/example.imzML", 1024, 2),
            # ("data/example.imzML", 1024, 8),
            # ("data/example.imzML", 1024, 32),
            # ("data/example.imzML", 1024, 64),
            # # batch_size = 2048
            # ("data/example.imzML", 2048, 1),
            # ("data/example.imzML", 2048, 2),
            # ("data/example.imzML", 2048, 8),
            # ("data/example.imzML", 2048, 32),
            # ("data/example.imzML", 2048, 64),
            # # batch_size = 4096
            # ("data/example.imzML", 4096, 1),
            # ("data/example.imzML", 4096, 2),
            # ("data/example.imzML", 4096, 8),
            # ("data/example.imzML", 4096, 32),
            # ("data/example.imzML", 4096, 64),
        ],
    )
    @pytest.mark.benchmark(timer=time.perf_counter)
    def test_batch_reading_benchmark(self, benchmark, data_file_path, batch_size, max_workers):
        logger.info("Starting data batch reading benchmark")
        benchmark.pedantic(
            data_batch_reading_full_benchmark,
            args=(data_file_path, batch_size, max_workers),
            rounds=ROUND,
            iterations=1,
            warmup_rounds=1,
        )
        drop_caches()
        logger.info("Finished data batch reading benchmark")
