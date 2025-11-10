from numpy import kaiser
from module.msi_module import MSI
from module.ms_module import MS
from module.msi_data_manager_msi import MSIDataManager  # 导入 .msi 文件的加载器
from module.msi_data_manager_zys import MSIDataManagerZYS
from preprocess.ms_preprocess import MSIPreprocessor
from logger import get_logger
logger = get_logger("example")

# Run examples when executing this file directly
if __name__ == "__main__":
    from module.ms_data_manager_imzml import MSDataManagerImzML
    from module.ms_module import MS
    FILE_PATH = "data/example.imzML"
    ms = MS()
    ms_md = MSDataManagerImzML(ms, filepath=FILE_PATH)
    ms_md.load_full_data_from_file()
    sp = ms[0]

    denoised = MSIPreprocessor.noise_reduction(
        data=sp,
        method="savgol",
        window=10,
        #window default = 5
        polyorder=3
        #polyorder default = 2
    )

    # Plotting
    denoised.plot(
        figsize=(12, 8),
        dpi=300,
        mz_range=(500.0, 510.0),
        intensity_range=(0.0, 1.5),
        original=sp,
        title_suffix='savgol'
    )
