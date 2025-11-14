from module.ms_module import MS
from preprocess.ms_preprocess import MSIPreprocessor
from module.ms_data_manager_imzml import MSDataManagerImzML
from tools.plot import plot_spectrum
from logger import get_logger

logger = get_logger("example")

# Run examples when executing this file directly
if __name__ == "__main__":
    FILE_PATH = "data/example.imzML"
    ms = MS()
    ms_md = MSDataManagerImzML(ms, filepath=FILE_PATH)
    ms_md.load_full_data_from_file()
    ms_md.inspect_data()

    for i, spectrum in enumerate(ms):
        denoised = MSIPreprocessor.noise_reduction_spectrum(spectrum,method="ma")
        ms[i].mz_list = denoised.mz_list
        ms[i].intensity = denoised.intensity

    denoise_test = ms[0]
    plot_spectrum(denoise_test,
                  mz_range=(500.0, 510.0),
                  intensity_range=(0, 1.2),)
    normalized_unit = MSIPreprocessor.normalization_spectrum(
    data=denoised,
    method="tic",          # or "median"
    scale_method="unit"    # min-max scaling to [0, 1]
    )

    plot_spectrum(
            base = normalized_unit,
            mz_range=(500.0, 510.0),
            intensity_range=(0.0, 0.1),
            title_suffix='TIC_normalized_unit'
        )