
from massflow.data_manager.ms_data_manager_imzml import MSDataManagerImzML
from massflow.module.mass_spectrum_set import MassSpectrumSet
from massflow.preprocess.dm_pre_fun import Preprocess
from massflow.tools.logger import get_logger
logger = get_logger("massflow")

def main():
    file_path = "data/example.imzML"

    mass_data = MassSpectrumSet()

    with MSDataManagerImzML(mass_data, filepath=file_path) as data_manager:
        data_manager.load_head_data()

        logger.info("Start async pipeline noise reduction with ma_loop")
        processed_data_manager = (
            Preprocess.pipeline(
                data_manager=data_manager,
                batch_size=128,
                temp_dir="./temp_noise_async_data",
                queue_ab_size=5,
                queue_bc_size=10,
                keep_order=False,
            )
            .noise_reduction(
                method="ma_loop",
                window=10,
                numba_max_threads=10,
            )
            .start()
        )

        logger.info(f"Async pipeline noise reduction finished. spectra={len(processed_data_manager.ms)}")
if __name__ == "__main__":
    main()
