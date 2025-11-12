from tkinter import X
from module.ms_module import MS
from preprocess.ms_preprocess import MSIPreprocessor
from module.ms_data_manager_imzml import MSDataManagerImzML
from tools.plot import plot_spectrum
from logger import get_logger

logger = get_logger("example")

if __name__ == "__main__":
    FILE_PATH = "data/example.imzML"
    ms = MS()
    ms_md = MSDataManagerImzML(ms, filepath=FILE_PATH)
    ms_md.load_full_data_from_file()
    ms_md.inspect_data()
    sp = ms[0]

    denoised = MSIPreprocessor.noise_reduction_spectrum(sp, method="ma")

    # denoised_2 = MSIPreprocessor.noise_reduction_spectrum(denoised,method="wavelet")
    # Plotting
    plot_spectrum(
        base=sp,
        target=denoised,
        mz_range=(400, 420),
        intensity_range=(-1, 10),
        metrics_box=True,
    )
