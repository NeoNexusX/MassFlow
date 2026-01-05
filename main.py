
from massflow.module.ms_data_manager_imzml import MSDataManagerImzML
from massflow.module.mass_spectrum_set import MassSpectrumSet
from massflow.logger import get_logger
from massflow.tools.plot import plot_spectrum
logger = get_logger(__name__)

def main():
    logger.info("Hello from massflow!")
    FILE_PATH = "data/example.imzML"
    ms = MassSpectrumSet()
    with MSDataManagerImzML(ms=ms, target_locs=[(1, 1), (50, 50)], filepath=FILE_PATH) as manager:
        manager.load_full_data_from_file()
        manager.inspect_data()
        spectrum = ms[0]

        # Use raw file descriptor for swapping to disk
        spectrum.swaper.swap_out2disk(manager.swap_fd) # type: ignore
        logger.info(f"Intensity: {spectrum.intensity}")
        input("Press Enter to continue...")


if __name__ == "__main__":
    main()