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
    # Create MS collection and manager
    # Auto-created from parser
    with MSDataManagerImzML(ms=ms,target_locs=[(1, 1), (50, 50)],filepath=FILE_PATH) as manager:

        # Load data with lazy-loading placeholders
        manager.load_full_data_from_file()
        spectrum = ms[0]
        spectrum.plot()
