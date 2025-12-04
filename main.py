
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
    
    spectrum_1 = ms[0]
    mz_range = (400, 500)
    spectrum_2 = ms[1]
    plot_spectrum(base = spectrum_1,
                  target = spectrum_2,
                  mz_range=mz_range)

if __name__ == "__main__":
    main()