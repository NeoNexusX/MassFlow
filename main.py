from massflow.data_manager import MSDataManagerImzML
from massflow.tools import get_logger
from massflow.module import MassSpectrumSet
from massflow.preprocess.flat_pre_fun import FlatPreprocess
from massflow.tools.plot import plot_spectrum

logger = get_logger("massflow")


def main():
    file_path = "data/example.imzML"

    with MSDataManagerImzML(filepath=file_path) as data_manager:
        data_manager.load_head_data()

        # Use flat_generator + flat numba normalization for faster compute path.
        processed_data_manager = MSDataManagerImzML(MassSpectrumSet(), temp_dir="temp")
        processed_data_manager.copy_meta(data_manager)

        for mz_flat, intensity_flat, lengths, coordinates in data_manager.flat_generator(
            batch_size=512,
            include_mz=True,
            max_threads=16,
        ):
            normalized_flat = FlatPreprocess.normalization_flat(
                intensity=intensity_flat,
                method="tic_numba",
                lengths=lengths,
            )
            processed_data_manager.swap_flat_data_out2disk(
                mz_flat=mz_flat,
                intensity_flat=normalized_flat,
                lengths=lengths,
                coordinates=coordinates,
            )

        processed_data_manager.close_writer()
        processed_data_manager.load_head_data()

        logger.info(
            "Normalization finished with flat_generator (method=tic_numba). "
            f"Processed spectra: {len(processed_data_manager.ms)}"
        )

        if len(processed_data_manager.ms) == 0:
            logger.warning("No spectra found in the processed dataset")
            processed_data_manager.close()
            return

        processed_spectrum = processed_data_manager.ms[0]
        plot_spectrum(
            base=processed_spectrum,
            mz_range=(200, 400),
            metrics_box=False,
        )

        processed_data_manager.close()

if __name__ == "__main__":
    main()
