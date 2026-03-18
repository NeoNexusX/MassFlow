import time
import os
from functools import partial
import pytest
from massflow.module.spectrum_imzml import SpectrumImzML
from massflow.data_manager.ms_data_manager_imzml import MSDataManagerImzML
from massflow.tools.logger import get_logger

logger = get_logger("test_swap")

ROUNDS = 5


def data_swap_out_direct(raw_dm, limit=-1):
    """
    Test: In-memory data -> Write to disk (Swap Out), using writer directly
    """

    limit = len(raw_dm.ms) if limit == -1 else limit

    swap_dm = MSDataManagerImzML(temp_dir="./temp_swap")
    swap_dm.copy_meta(raw_dm)
    writer = swap_dm.writer

    for i, raw_spec in enumerate(raw_dm.ms):
        if i >= limit:
            break

        writer.add_spectrum(
            raw_spec.mz_list, raw_spec.intensity, raw_spec.coordinate.get_tuple()
        )

    swap_dm.close_writer()
    return swap_dm


def data_swap_out_object(raw_dm, limit=-1):
    """
    Test: In-memory data -> Write to disk (Swap Out), using SpectrumImzML object method
    """

    limit = len(raw_dm.ms) if (limit is None or limit == -1) else limit

    dm_swap = MSDataManagerImzML(temp_dir="./temp_swap")
    dm_swap.copy_meta(raw_dm)
    writer = dm_swap.writer

    for i, old_spec in enumerate(raw_dm.ms):
        if i >= limit:
            break

        spec_processed_imzml = SpectrumImzML(
            index=i,
            mz_list=old_spec.mz_list,
            intensity=old_spec.intensity,
            coordinates=old_spec.coordinate,
        )

        spec_processed_imzml.swap_out2disk(writer)

    dm_swap.close_writer()
    return dm_swap


def data_swap_in_benchmark(swap_dm, mode="multithread"):
    """
    Test: Disk file -> Load and access in memory (Swap In)
    """
    dm = swap_dm
    dm.load_full_data_from_file()
    dm.inspect_data()

    if mode == "multithread":
        batch_generator = dm.get_batch_generator(batch_size=512, max_threads=32)

        for batch in batch_generator:
            assert len(batch[1]) > 0

            # clean batch memory
            for spec in batch:
                spec.clear_data()
    else:

        # read in memory
        for spec in dm.ms:
            assert spec.intensity is not None

        # clean memory
        for spec in dm.ms:
            spec.clear_data()

    logger.info("Completed loading and accessing swapped data.")


class TestDataSwap:
    """
    Disk swap performance test class
    """

    @pytest.fixture(scope="class")
    def ms_raw_data(self):

        filepath = r"./data/example.imzML"

        if not os.path.exists(filepath):
            pytest.skip(f"Original data not found at {filepath}")

        # load meta data
        dm_raw = MSDataManagerImzML(filepath=filepath)
        dm_raw.load_head_data()
        dm_raw.inspect_data()
        # load data check
        batch_generator = dm_raw.batch_generator()
        for batch in batch_generator:
            assert len(batch[1]) > 0

        return dm_raw

    @pytest.mark.benchmark(timer=time.perf_counter)
    def test_swap_out_performance_direct(self, benchmark, ms_raw_data, limit=-1):
        """
        Performance test: Data write speed to disk
        """
        raw_dm = ms_raw_data

        logger.info(f"Starting Swap Out benchmark for {limit} spectra")

        dm_result = benchmark.pedantic(
            data_swap_out_direct,
            args=(raw_dm, limit),
            rounds=ROUNDS,
            iterations=1,
            warmup_rounds=1,
        )

        dm_result.close()
        logger.info("Finished Swap Out benchmark")

    @pytest.mark.benchmark(timer=time.perf_counter)
    def test_swap_out_performance_object(self, benchmark, ms_raw_data, limit=-1):
        """
        Performance test: Data write speed to disk
        """
        raw_dm = ms_raw_data

        logger.info(f"Starting Swap Out benchmark for {limit} spectra")

        dm_result = benchmark.pedantic(
            data_swap_out_object,
            args=(raw_dm, limit),
            rounds=ROUNDS,
            iterations=1,
            warmup_rounds=1,
        )

        dm_result.close()
        logger.info("Finished Swap Out benchmark")

    @pytest.mark.benchmark(timer=time.perf_counter)
    @pytest.mark.parametrize("mode", ["normal", "multithread"])
    def test_swap_in_performance(self, benchmark, ms_raw_data, mode, limit=-1):
        """
        Performance test: Data read speed from disk
        """
        logger.info(f"Starting Swap In benchmark for {limit} spectra")

        raw_dm = ms_raw_data

        def set_up_swap_data():
            dm_swap_new = data_swap_out_direct(raw_dm, limit=limit)
            raw_dm.clear_all_data_memory()
            return (dm_swap_new,), {}

        benchmark.pedantic(
            partial(data_swap_in_benchmark, mode=mode),
            setup=set_up_swap_data,
            rounds=ROUNDS,
            iterations=1,
            warmup_rounds=1,
        )

        logger.info("Finished Swap In benchmark")
