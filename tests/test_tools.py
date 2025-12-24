import time
import functools
from massflow.logger import get_logger
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
