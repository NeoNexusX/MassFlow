
from massflow.data_manager.ms_data_manager_imzml import MSDataManagerImzML
from massflow.module.mass_spectrum_set import MassSpectrumSet
from massflow.tools.logger import get_logger
logger = get_logger(__name__)

def main():
    logger.info("Hello from massflow!")
    file_path = "data/example.imzML"
    ms = MassSpectrumSet()
    with MSDataManagerImzML(ms=ms, target_locs=[(1, 1), (50, 50)], filepath=file_path) as manager:
        manager.load_head_data()
        manager.inspect_data()
        spectrum = ms[0]

        logger.info(f"Intensity: {spectrum.intensity}")
        input("Press Enter to continue...")


if __name__ == "__main__":
    main()
