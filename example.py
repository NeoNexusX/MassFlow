from module.ms_module import MS
from preprocess.ms_preprocess import MSIPreprocessor
from module.ms_data_manager_imzml import MSDataManagerImzML
from tools.plot import plot_spectrum
from logger import get_logger

logger = get_logger("example")

if __name__ == "__main__":
    FILE_PATH = "data/example2.imzML"
    ms = MS()
    ms_md = MSDataManagerImzML(ms, filepath=FILE_PATH)
    ms_md.load_full_data_from_file()
    ms_md.inspect_data()
    sp = ms[100]

    corrected_sp, baseline = MSIPreprocessor.baseline_correction(
        data=sp,
        method="asls",
        lam=5e7,
        p=0.008,
        niter=25,
        baseline_scale=0.9
    )

    # Plotting
    plot_spectrum(
    base=sp,
    target=corrected_sp,
    mz_range=(40, 100),
    intensity_range=(-1, 100),
    metrics_box=True,
    title_suffix="ASLS"
    )