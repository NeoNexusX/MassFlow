import time
import pytest
import os
from massflow.module.spectrum_imzml import SpectrumImzML
from massflow.module.mass_spectrum_set import MassSpectrumSet
from massflow.module.ms_data_manager_imzml import MSDataManagerImzML
from massflow.logger import get_logger

logger = get_logger("test_swap")

ROUND = 1

def ms_raw_data():
    """
    create data
    """
    filepath = r"./data/example.imzML"

    if not os.path.exists(filepath):
        pytest.skip(f"Original data not found at {filepath}")
        
    ms_raw = MassSpectrumSet()
    dm_raw = MSDataManagerImzML(ms_raw, filepath=filepath)

    return ms_raw, dm_raw

def ms_raw_data_readin(dm_raw):
    """
    Fixture: Load raw data, shared for all test cases
    """
    dm_raw.load_full_data_from_file()
    batch_generator = dm_raw.get_batch_generator()
    for batch in batch_generator:
        pass

    return dm_raw

def data_swap_out_direct(ms_raw, limit=-1):
    """
    Test: In-memory data -> Write to disk (Swap Out), using writer directly
    """

    limit = len(ms_raw) if limit == -1 else limit

    ms_swap = MassSpectrumSet()
    dm_swap = MSDataManagerImzML(ms_swap,temp_dir="./temp_swap")
    dm_swap.copy_meta(ms_raw)
    writer = dm_swap.writer

    for i, spec_raw in enumerate(ms_raw):
        if i >= limit:
            break

        writer.addSpectrum(spec_raw.mz_list, 
                           spec_raw.intensity, 
                           spec_raw.coordinates.get_tuple())
        #clean old memory
        spec_raw.clear_data()

    dm_swap.close_writer()
    return dm_swap
    
def data_swap_out_object(ms_raw, limit=-1):
    """
    Test: In-memory data -> Write to disk (Swap Out), using SpectrumImzML object method
    """

    limit = len(ms_raw) if (limit is None or limit == -1) else limit

    ms_swap = MassSpectrumSet()
    dm_swap = MSDataManagerImzML(ms_swap,temp_dir="./temp_swap")
    dm_swap.copy_meta(ms_raw)
    writer = dm_swap.writer

    for i, old_spec in enumerate(ms_raw):
        if i >= limit:
            break

        spec_processed_imzml = SpectrumImzML(index=i, 
                                             mz_list=old_spec.mz_list, 
                                             intensity=old_spec.intensity, 
                                             coordinates=old_spec.coordinates)

        spec_processed_imzml.swap_out2disk(writer)

        #clean raw data:
        old_spec.clear_data()

    dm_swap.close_writer()
    return dm_swap

def data_swap_in_benchmark(dm_swap):
    """
    Test: Disk file -> Load and access in memory (Swap In)
    """
    dm = dm_swap
    dm.load_full_data_from_file()
    batch_generator = dm.get_batch_generator(batch_size=256, max_threads=8)

    for batch in batch_generator:
        assert len(batch[1]) >0
    logger.info("Completed loading and accessing swapped data.")

    # input("stage 4: load swap data")

    # Cleanup to release memory
    for spec in dm.ms:
        spec.clear_data()

    # input("stage 5: clean swap data")

def set_up_swap_data():
    ms_raw, dm_raw = ms_raw_data()  # Unpack the tuple
    # input("stage 1: load raw meta data")
    ms_raw_data_readin(dm_raw) # load data into ms_raw
    # input("stage 2: load raw data")
    dm_swap_new = data_swap_out_direct(ms_raw, limit=-1)
    # input("stage 3: write swap data to disk, raw data + ")

    return (dm_swap_new, ), {}

def pipeline_swap_data():
    pass

class TestDataSwap:
    """
    Disk swap performance test class
    """

    @pytest.mark.benchmark(timer=time.perf_counter)
    def test_swap_out_performance_direct(self, benchmark, ms_raw_data, limit=-1):
        """
        Performance test: Data write speed to disk
        """                 
        ms_raw, _ = ms_raw_data

        logger.info(f"Starting Swap Out benchmark for {limit} spectra")

        dm_result = benchmark.pedantic(data_swap_out_direct,
                                       args=(ms_raw, limit),
                                       rounds=ROUND,
                                       iterations=1,
                                       warmup_rounds=1)

        dm_result.close()
        logger.info("Finished Swap Out benchmark")

    @pytest.mark.benchmark(timer=time.perf_counter)
    def test_swap_out_performance_object(self, benchmark, ms_raw_data, limit=-1):
        """
        Performance test: Data write speed to disk
        """                 
        ms_raw, _ = ms_raw_data

        logger.info(f"Starting Swap Out benchmark for {limit} spectra")

        dm_result = benchmark.pedantic(data_swap_out_object,
                                       args=(ms_raw, limit),
                                       rounds=ROUND,
                                       iterations=1,
                                       warmup_rounds=1)

        dm_result.close()
        logger.info("Finished Swap Out benchmark")

    @pytest.mark.benchmark(timer=time.perf_counter)
    def test_swap_in_performance(self, benchmark, limit=-1):
        """
        Performance test: Data read speed from disk
        """

        logger.info(f"Starting Swap In benchmark for {limit} spectra")

        benchmark.pedantic(data_swap_in_benchmark,
                           setup=set_up_swap_data,
                           rounds=ROUND,
                           iterations=1,
                           warmup_rounds=0)

        logger.info("Finished Swap In benchmark")
