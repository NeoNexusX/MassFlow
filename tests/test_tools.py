import time
import functools
import pytest
from massflow.data_manager import MSDataManagerImzML
from massflow.tools.logger import get_logger
from massflow.tools.infer_spectrum_type import resolve_spectrum_type, infer_spectrum_type

logger = get_logger("test_tools")

#暂时没用到
def time_calculator(test_times=1):

    def timer(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            times = []
            for i in range(test_times):
                t0 = time.perf_counter()
                func(*args, **kwargs)
                t1 = time.perf_counter()
                times.append(t1 - t0)
                logger.info(f"test times : {i}, time : {t1 - t0:.6f} s")

            logger.info(f"total time : {sum(times):.6f} s ,test times : {test_times})"
                        f"average time : {sum(times)/test_times:.6f} s")
        return wrapper

    return timer

@pytest.mark.parametrize("data_file_path", ["data/example_file.imzML"])
def test_infer_spectrum_type(data_file_path: str) -> None:
    """test infer_spectrum_type function with example imzML file"""
    raw_dm = MSDataManagerImzML(filepath=data_file_path)
    raw_dm.load_head_data()

    raw_spectrum_type = resolve_spectrum_type(raw_dm)
    inferred_spectrum_type = infer_spectrum_type(raw_dm)

    try:
        if raw_spectrum_type != inferred_spectrum_type:
            pytest.fail(
                f"Spectrum type inference failed for {data_file_path}: Expected: {raw_spectrum_type}, Inferred: {inferred_spectrum_type}"
            )
        else:
            logger.info(
                f"Spectrum type inference passed for {data_file_path}: Expected: {raw_spectrum_type}, Inferred: {inferred_spectrum_type}"
            )
    finally:
        raw_dm.close()
