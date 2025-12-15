
from massflow.module.ms_data_manager_imzml import MSDataManagerImzML
from massflow.module.ms_module import MS
from massflow.logger import get_logger
from massflow.tools.plot import plot_spectrum
logger = get_logger(__name__)

def main():
    logger.info("Hello from massflow!")
    FILE_PATH = "data/example.imzML"
    ms = MS()
    with MSDataManagerImzML(ms=ms, target_locs=[(1, 1), (50, 50)], filepath=FILE_PATH) as manager:
        manager.load_full_data_from_file()
        manager.inspect_data()
        


if __name__ == "__main__":
    main()