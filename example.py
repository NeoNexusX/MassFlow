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
    sp = ms[0]

    denoised = MSIPreprocessor.noise_reduction(
        data=sp,
        method="savgol",
        window=5,
        polyorder=2
    )

    # Plotting
    plot_spectrum(
        base=sp,
        target=denoised,
        mz_range=(500.0, 505.0),
        intensity_range=(0, 1.2),
        title_suffix='savgol',
        overlay=False,
    )