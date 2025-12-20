import time
import pytest
import os
from massflow.module.spectrum_imzml import SpectrumImzML
from massflow.module.mass_spectrum_set import MassSpectrumSet
from massflow.module.ms_data_manager_imzml import MSDataManagerImzML
from massflow.logger import get_logger

logger = get_logger("test_swap")

ROUND = 10

def data_swap_out_direct(ms_raw, limit=-1):
    """
    测试：内存数据 -> 磁盘写入 (Swap Out), 直接使用 writer
    """

    limit = len(ms_raw) if limit == -1 else limit

    ms_swap = MassSpectrumSet()
    dm_swap = MSDataManagerImzML(ms_swap)
    dm_swap.copy_meta(ms_raw)
    writer = dm_swap.writer
    
    for i, spec in enumerate(ms_raw):
        if i >= limit:
            break
        spec_processed = spec

        writer.addSpectrum(spec_processed.mz_list, 
                           spec_processed.intensity, 
                           spec_processed.coordinates.get_tuple())

    dm_swap.close_writer()
    return dm_swap

def data_swap_out_object(ms_raw, limit=-1):
    """
    测试：内存数据 -> 磁盘写入 (Swap Out), 使用 SpectrumImzML 对象方法
    """
    ms_swap = MassSpectrumSet()
    dm_swap = MSDataManagerImzML(ms_swap,temp_dir="./temp_swap")
    dm_swap.copy_meta(ms_raw)
    writer = dm_swap.writer

    limit = len(ms_raw) if (limit is None or limit == -1) else limit

    for i, spec in enumerate(ms_raw):
        if i >= limit:
            break
        spec_processed = spec

        spec_processed_imzml = SpectrumImzML(
            index=i, 
            mz_list=spec_processed.mz_list, 
            intensity=spec_processed.intensity, 
            coordinates=spec.coordinates
        )

        spec_processed_imzml.swap_out2disk(writer)

    dm_swap.close_writer()
    return dm_swap

def data_swap_in_benchmark(dm_swap):
    """
    测试：磁盘文件 -> 内存加载与访问 (Swap In)
    """
    dm_swap.load_full_data_from_file()
    
    batch_generator = dm_swap.get_batch_generator()
    for batch in batch_generator:
        for spec in batch:
            _ = spec.mz_list
            _ = spec.intensity

class TestDataSwap:
    """
    磁盘交换性能测试类
    """
    
    @pytest.fixture(scope="class")
    def ms_raw_data(self):
        """
        固件：加载原始数据，供所有测试用例共享
        """
        filepath = r"./data/example.imzml"

        if not os.path.exists(filepath):
            pytest.skip(f"Original data not found at {filepath}")
            
        ms_raw = MassSpectrumSet()
        dm_raw = MSDataManagerImzML(ms_raw, filepath=filepath)
        dm_raw.load_full_data_from_file()
        batch_generator = dm_raw.get_batch_generator()
        for batch in batch_generator:
            for spec in batch:
                # 强制触发懒加载并验证数据有效
                _ = spec.mz_list
                _ = spec.intensity

        return ms_raw, dm_raw
    

    @pytest.mark.parametrize("limit", [-1])
    @pytest.mark.parametrize("method",['object','direct'])
    @pytest.mark.benchmark(timer=time.perf_counter)
    def test_swap_out_performance(self, benchmark, ms_raw_data, limit, method):
        """
        性能测试：数据写入磁盘的速率
        """
        out_func = data_swap_out_direct if method == "direct" else data_swap_out_object
        ms_raw, dm_raw = ms_raw_data

        logger.info(f"Starting Swap Out benchmark for {limit} spectra")

        dm_result = benchmark.pedantic(
            out_func,
            args=(ms_raw, limit),
            rounds=ROUND,
            iterations=1,
            warmup_rounds=1
        )

        dm_result.close()
        logger.info("Finished Swap Out benchmark")

    @pytest.mark.parametrize("limit", [-1])
    @pytest.mark.benchmark(timer=time.perf_counter)
    def test_swap_in_performance(self, benchmark, ms_raw_data, limit):
        """
        性能测试：数据从磁盘回读的速率
        """
        ms_raw, dm_raw = ms_raw_data  # Unpack the tuple

        def set_up_swap_data():
            dm_swap_new = data_swap_out_direct(ms_raw, limit=limit)
            return (dm_swap_new, ), {}

        logger.info(f"Starting Swap In benchmark for {limit} spectra")

        benchmark.pedantic(
            data_swap_in_benchmark,
            setup=set_up_swap_data,
            rounds=ROUND,
            iterations=1,
            warmup_rounds=1
        )

        logger.info("Finished Swap In benchmark")