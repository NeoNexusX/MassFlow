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
    FILE_PATH = "data/20250421_20um.imzML"
    #example for read imzml data as ms ：
    ms = MS()
    ms_dm = MSDataManagerImzML(ms, filepath=FILE_PATH)
    ms_dm.load_full_data_from_file()
    ms_dm.inspect_data()
    spectrum1 = ms[0]
    spectrum2 = MSIPreprocessor.noise_reduction(spectrum1, method='ma_ns', window=10)

    # example usage for MSIDataManager:
    # msi = MSI(name='example1', version=1.0, mask=None, need_base_mask=True)
    # msi_dm = MSIDataManagerMSI(msi,filepath=FILE_PATH)
    # msi_dm.load_full_data_from_file()
    # msi_dm.inspect_data()
    # assert msi.meta.mz_num >0 , "meta_mz_num must be greater than 0"

    # # example for MSI_DataManager read file with specific mz range file:
    # msi2 = MSI(name='example2', version=1.0, need_base_mask=True)
    # msi_dm2 = MSIDataManagerMSI(msi2, target_mz_range=[100, 150.5], filepath=FILE_PATH)
    # msi_dm2.inspect_data()
    # msi_dm2.load_full_data_from_file()
    # msi2.plot_msi()


    # example for MSI_DataManager_ZYS read file with specific mz range file:
    # msi3 = MSI(name='20250421_6L_1uL_0,2MPa_3000V_95ACN_10um01.mat', version=1.0)
    # FILE_PATH_ZYS = "./data/20250421_6L_1uL_0,2MPa_3000V_95ACN_10um01.mat"
    # msi_dm_zys = MSIDataManagerZYS(msi3, target_mz_range=[100, 170.5], filepath = FILE_PATH_ZYS)
    # msi_dm_zys.load_full_data_from_file()
    # msi_dm_zys.rebuild_hdf5_file_from_zys()
    # msi_dm_zys.inspect_data()
    # # # example usage for MSI_DataManager_ZYS to write a new hdf5 file:
    # msi_dm_zys.write2local(mode='merge',output_fold='./data')
    # msi_dm_zys.write2local(mode='split')

    # # example usage for get one msi image by mz value:
    # msi_list = msi3.get_msi_by_mz(mz_value_min=88.1122, tol=1e-3)
    # image = msi_list[0].msroi.T
    # base_mask  = msi_list[0].base_mask.T
