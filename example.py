from module.msi_module import MSI
from module.ms_module import MS
from module.msi_data_manager_msi import MSIDataManagerMSI  # 导入 .msi 文件的加载器
from module.msi_data_manager_zys import MSIDataManagerZYS
from module.ms_data_manager_imzml import MSDataManagerImzML
from preprocess.ms_preprocess import MSIPreprocessor
from logger import get_logger
logger = get_logger("example")

# Run examples when executing this file directly
if __name__ == "__main__":

    # Replace with your data file path
    FILE_PATH = "data/example.imzML"
    #example for read imzml data as ms ：
    ms = MS()
    ms_dm = MSDataManagerImzML(ms, filepath=FILE_PATH)
    ms_dm.load_full_data_from_file()
    ms_dm.inspect_data()
    spectrum1 = ms[0]
    print(min(spectrum1.intensity))
    spectrum1.plot()
    spectrum2 = MSIPreprocessor.noise_reduction(spectrum1, method='wavelet')
    spectrum2.plot()
