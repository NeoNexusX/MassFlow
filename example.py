from module.msi_module import MSI
from module.ms_module import MS
from module.msi_data_manager_msi import MSIDataManager  # 导入 .msi 文件的加载器
from module.msi_data_manager_zys import MSIDataManagerZYS
from module.ms_data_manager_imzml import MSDataManagerImzML
from preprocess.ms_preprocess import MSIPreprocessor
from logger import get_logger
logger = get_logger("example")

# Run examples when executing this file directly
if __name__ == "__main__":

    # Replace with your data file path
    # FILE_PATH = "data/MSI_default_merge_1.0.msi"
    # msi = MSI()
    # msi_md = MSIDataManager(msi, filepath=FILE_PATH)
    # msi_md.load_full_data_from_file()
    # msi_md.inspect_data()

    # # # Replace with your data file path
    # FILE_PATH_ZYS = "data/2Moderate_20250325_15uL_90ACN+HAC_4000V_20um.mat"
    # msi_zys = MSI()
    # msi_zys_md = MSIDataManagerZYS(msi_zys, filepath=FILE_PATH_ZYS)
    # msi_zys_md.load_full_data_from_file()
    # msi_zys_md.inspect_data()
    # msi_zys_md.write2local(output_fold="./data")

    FILE_PATH = "data/example.imzML"
    ms = MS()
    ms_md = MSDataManagerImzML(ms, filepath=FILE_PATH)
    ms_md.load_full_data_from_file()
    ms_md.inspect_data()
    ms.plot_ms_mask()
    spectrum = ms[0]
    spectrum.plot()