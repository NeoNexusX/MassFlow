from massflow.data_manager import MSDataManagerImzML
from massflow.module import MassSpectrumSet
from massflow.preprocess import Preprocessor
from massflow.tools import get_logger
logger = get_logger("massflow")

def main():
    file_path = "data/example.imzML"

    mass_data = MassSpectrumSet()

    with MSDataManagerImzML(mass_data, filepath=file_path) as data_manager:
        data_manager.load_head_data()

        logger.info("Start async pipeline noise reduction with ma_loop")
        processed_data_manager = (
            Preprocessor(data_manager, numba_max_threads=10)
            .noise_reduction(
                method="ma_loop",
                window=10,
            )
            .start()
        )

        logger.info(f"Async pipeline noise reduction finished. spectra={len(processed_data_manager.ms)}")
if __name__ == "__main__":
    main()
