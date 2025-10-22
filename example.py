from module.msi_data_manager_msi import MSIDataManagerMSI
from module.msi_data_manager_zys import MSIDataManagerZYS
import os
from module.msi_module import MSI
from module.msi_data_manager_imzml import MSIDataManagerImzML

# Run examples when executing this file directly
if __name__ == "__main__":


    # Replace with your data file path
    FILE_PATH = "./data/MSI_2Moderate_20250325_15uL_90ACN+HAC_4000V_20um_merge_1.0.msi"
    # example usage for MSIDataManager:
    msi = MSI(name='example1', version=1.0, mask=None, need_base_mask=True)
    msi_dm = MSIDataManagerMSI(msi,filepath=FILE_PATH)
    msi_dm.load_full_data_from_file()
    msi_dm.inspect_data()
    assert msi.meta_mz_num >0 , "meta_mz_num must be greater than 0"

    # example for MSI_DataManager read file with specific mz range file:
    msi2 = MSI(name='example2', version=1.0, mask=None, need_base_mask=True)
    msi_dm2 = MSIDataManagerMSI(msi2, target_mz_range=[101, 102], filepath=FILE_PATH)
    msi_dm2.inspect_data()
    msi_dm2.load_full_data_from_file()
    msi2.plot_msi()

    # example for MSI_DataManager_ZYS read file with specific mz range file:
    msi3 = MSI(name='2Moderate_20250325_15uL_90ACN+HAC_4000V_20um', version=1.0, mask=None, need_base_mask=True)
    FILE_PATH_ZYS = "./data/2Moderate_20250325_15uL_90ACN+HAC_4000V_20um.mat"
    msi_dm_zys = MSIDataManagerZYS(msi3, filepath = FILE_PATH_ZYS)
    msi_dm_zys.load_full_data_from_file()
    msi_dm_zys.rebuild_hdf5_file_from_zys()
    msi_dm_zys.inspect_data()
    # example usage for MSI_DataManager_ZYS to write a new hdf5 file:
    msi_dm_zys.write2local(mode='merge',output_fold='./data')
    msi_dm_zys.write2local(mode='split')

    # example usage for get one msi image by mz value:
    msi_list = msi3.get_msi_by_mz(mz_value_min=88.1122, tol=1e-3)
    image = msi_list[0].msroi.T
    base_mask  = msi_list[0].base_mask.T
