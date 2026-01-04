import time
import asyncio
import pytest
from massflow.module.ms_data_manager_imzml import MSDataManagerImzML
from massflow.logger import get_logger

logger = get_logger("test_read")

ROUND = 1


# -----------------------------------------------------------------------------
# Async Benchmark Helper
# -----------------------------------------------------------------------------
async def _async_runner(data_file_path, batch_size, max_concurrency):
    """
    Benchmark for asyncio-based get_batch_generator_async method.
    """
    # Use context manager for automatic resource cleanup (ms_md.close())
    with MSDataManagerImzML(filepath=data_file_path) as ms_md:
        # 1. Load metadata and placeholders (synchronous method)
        ms_md.load_full_data_from_file()

        # 2. Get the async generator
        batch_generator = ms_md.get_batch_generator_async(
            batch_size=batch_size, max_concurrency=max_concurrency
        )

        # 3. Asynchronous iteration (core part of the test)
        async for batch in batch_generator:
            # Simple verification that data has been loaded
            if len(batch) > 0:
                # Verify intensity has been preloaded (should not be None)
                assert len(batch[0].intensity) >= 0


def data_reading_benchmark(data_file_path="data/example.imzML"):
    with MSDataManagerImzML(filepath=data_file_path) as ms_md:
        ms_md.load_full_data_from_file()
        for spectrum_data in ms_md.ms:
            assert spectrum_data.intensity is not None


def data_batch_reading_async_benchmark(
    data_file_path="data/example.imzML", batch_size=1024, max_concurrency=2
):
    """
    Synchronous wrapper for benchmarking.
    pytest-benchmark works best with synchronous functions.
    We use asyncio.run() here to measure the full execution time of the async loop.
    """
    asyncio.run(_async_runner(data_file_path, batch_size, max_concurrency))


def data_batch_reading_full_benchmark(
    data_file_path="data/example.imzML", batch_size=1024, max_workers=2
):
    with MSDataManagerImzML(filepath=data_file_path) as ms_md:
        ms_md.load_full_data_from_file()
        ms_md.inspect_data()
        batch_generator = ms_md.get_batch_generator(
            batch_size=batch_size, max_threads=max_workers
        )
        for batch in batch_generator:
            assert len(batch) > 0


def data_batch_reading_benchmark_new(
    data_file_path="data/example.imzML", batch_size=1024, max_workers=2
):
    """
    Benchmark for the synchronous/threaded version (comparison baseline).
    """
    with MSDataManagerImzML(filepath=data_file_path) as ms_md:
        ms_md.load_full_data_from_file()
        ms_md.inspect_data()

        batch_generator = ms_md.get_batch_generator_new(
            batch_size=batch_size, max_threads=max_workers
        )
        for batch in batch_generator:
            if len(batch) > 0:
                assert batch[0].intensity is not None
                assert len(batch[0].intensity) >= 0


class TestDataReading:

    @pytest.mark.parametrize("data_file_path", ["data/example.imzML"])
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
        "data_file_path,batch_size,max_concurrency",
        [
            ("data/example.imzML", 256, 64),
            ("data/example.imzML", 256, 128),
            ("data/example.imzML", 1024, 128),
            ("data/example.imzML", 1024, 256),
            ("data/example.imzML", 2048, 256),
            ("data/example.imzML", 2048, 1024),
            ("data/example.imzML", 4096, 32),
            ("data/example.imzML", 4096, 64),
            ("data/example.imzML", 4096, 256),
            ("data/example.imzML", 4096, 512),
            ("data/example.imzML", 4096, 1024),
        ],
    )
    def test_batch_data_reading_benchmark_async(
        self, benchmark, data_file_path, batch_size, max_concurrency
    ):
        logger.info("Starting batch data reading async benchmark")

        # Use pedantic mode to run the synchronous wrapper
        benchmark.pedantic(
            data_batch_reading_async_benchmark,
            kwargs={
                "data_file_path": data_file_path,
                "batch_size": batch_size,
                "max_concurrency": max_concurrency,
            },
            rounds=ROUND,
            iterations=1,
            warmup_rounds=1,
        )
        logger.info("Finished batch data reading async benchmark")

    @pytest.mark.parametrize(
        "data_file_path,batch_size,max_workers",
        [
            ("data/example.imzML", 32, 2),
            ("data/example.imzML", 64, 32),
            ("data/example.imzML", 64, 64),
            ("data/example.imzML", 128, 8),
            ("data/example.imzML", 256, 32),
            ("data/example.imzML", 256, 64),
            ("data/example.imzML", 1024, 8),
            ("data/example.imzML", 4096, 32),
            ("data/example.imzML", 4096, 64),
        ],
    )
    def test_batch_data_reading_benchmark_new(
        self, benchmark, data_file_path, batch_size, max_workers
    ):
        logger.info("Starting batch data reading benchmark (new)")
        benchmark.pedantic(
            data_batch_reading_benchmark_new,
            args=(data_file_path, batch_size, max_workers),
            rounds=ROUND,
            iterations=1,
            warmup_rounds=1,
        )
        logger.info("Finished batch data reading benchmark (new)")

    @pytest.mark.parametrize(
        "data_file_path,batch_size,max_workers",
        [
            ("data/example.imzML", 32, 2),
            ("data/example.imzML", 64, 32),
            ("data/example.imzML", 64, 64),
            ("data/example.imzML", 128, 8),
            ("data/example.imzML", 256, 32),
            ("data/example.imzML", 256, 64),
            ("data/example.imzML", 1024, 8),
            ("data/example.imzML", 4096, 32),
            ("data/example.imzML", 4096, 64),
        ],
    )
    def test_batch_data_reading_benchmark_full(
        self, benchmark, data_file_path, batch_size, max_workers
    ):
        logger.info("Starting batch data reading full benchmark (original)")
        benchmark.pedantic(
            data_batch_reading_full_benchmark,
            args=(data_file_path, batch_size, max_workers),
            rounds=ROUND,
            iterations=1,
            warmup_rounds=1,
        )
