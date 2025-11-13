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

    smoother = MSIPreprocessor.noise_reduction_spectrum(sp,method="bi_ns")
    pickpeak = MSIPreprocessor.peak_pick_spectrum(smoother,return_type="area")

    # denoised_2 = MSIPreprocessor.noise_reduction_spectrum(denoised,method="wavelet")
    # Plotting
    plot_spectrum(
        base=sp,
        target=pickpeak,
        mz_range=(400, 420),
        intensity_range=(-1, 15),
        metrics_box=False,
        overlay=True,
        plot_mode=["line","stem"],
    )
