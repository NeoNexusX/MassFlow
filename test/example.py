from massflow.module.ms_module import MS
from massflow.preprocess.ms_preprocess import MSIPreprocessor
from massflow.module.ms_data_manager_imzml import MSDataManagerImzML
import time
from massflow.logger import get_logger
from typing import Optional
logger = get_logger("example")


def benchmark(test_times=1,
              mass_data: Optional[MS]=None,
              func=None):

    if mass_data is None or func is None:
        logger.error("mass_data or func is None, please check!")
        return

    #pre load data
    for spectrum_data in mass_data:
        _ = spectrum_data.mz_list

    # time start:
    times = []
    for i in range(test_times):
        t0 = time.perf_counter()
        # run function:
        func(mass_data)
        t1 = time.perf_counter()
        times.append(t1 - t0)
        logger.info(f"test times : {i}, time : {t1 - t0:.6f} s")

    logger.info(f"total time : {sum(times):.6f} s ,test times : {test_times})"
                f"average time : {sum(times)/test_times:.6f} s")


def denoise_benchmark(mass_spectrum: Optional[MS]=None):

    if mass_spectrum is not None:
        for i, spectrum in enumerate(mass_spectrum):
            denoised = MSIPreprocessor.noise_reduction_spectrum(spectrum,method="ma")
            mass_spectrum[i].mz_list = denoised.mz_list
            mass_spectrum[i].intensity = denoised.intensity



# Run examples when executing this file directly
if __name__ == "__main__":
    FILE_PATH = "data/example.imzML"
    ms = MS()
    with MSDataManagerImzML(ms, filepath=FILE_PATH) as ms_md:
        ms_md.load_full_data_from_file()
        ms_md.inspect_data()

        benchmark(test_times=1000,
                    mass_data=ms,
                    func=denoise_benchmark)