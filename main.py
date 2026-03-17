
from massflow.data_manager.ms_data_manager_imzml import MSDataManagerImzML
from massflow.module.mass_spectrum_set import MassSpectrumSet
from massflow.tools.plot import plot_spectrum
from massflow.tools.logger import get_logger
logger = get_logger("massflow")

def main():
    file_path = "data/15-0879 wb xenograft trp pathway dosed- rms_corrected.imzML"
    ms = MassSpectrumSet()
    with MSDataManagerImzML(ms=ms, target_locs=[(1, 1), (50, 50)], filepath=file_path) as manager:
        manager.load_head_data()
        spectrum = ms[0]

        plot_spectrum(spectrum,plot_mode=["stem"])

        logger.info(f"Intensity: {spectrum.intensity}")
        input("Press Enter to continue...")


if __name__ == "__main__":
    main()
