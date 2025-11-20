import time
from massflow.module.ms_module import MS
from massflow.logger import get_logger
logger = get_logger("__test_tools__")


def time_calculator(test_times=100):
    
    def timer(func):

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

# def denoise_benchmark(mass_spectrum: Optional[MS]=None):

#     if mass_spectrum is not None:
#         for i, spectrum in enumerate(mass_spectrum):
#             denoised = MSIPreprocessor.noise_reduction_spectrum(spectrum,method="ma")
#             mass_spectrum[i].mz_list = denoised.mz_list
#             mass_spectrum[i].intensity = denoised.intensity

# # Run examples when executing this file directly
# if __name__ == "__main__":
#     FILE_PATH = "data/example.imzML"
#     ms = MS()
#     with MSDataManagerImzML(ms, filepath=FILE_PATH) as ms_md:
#         ms_md.load_full_data_from_file()
#         ms_md.inspect_data()
