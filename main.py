from massflow.data_manager import MSDataManagerImzML
from massflow.tools import get_logger
from massflow.tools.plot import plot_spectrum
from massflow.tools.snr import calculate_snr_details, log_snr_details
from massflow.module import MassSpectrumSet
from massflow.preprocess.flat_pre_fun import FlatPreprocess
from pathlib import Path

logger = get_logger("massflow")


def plot_before_after(
    before_spectrum,
    after_spectrum,
    save_path: str | None = None,
) -> None:
    """Display overlay plot by default; save when a path is provided."""
    output_path = None
    if save_path:
        output_path = Path(save_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

    plot_spectrum(
        base=before_spectrum,
        target=after_spectrum,
        save_path=str(output_path) if output_path else None,
        overlay=True,
        metrics_box=True,
    )
    if output_path:
        logger.info(f"Saved before/after spectrum plot: {output_path}")
    else:
        logger.info("Displayed before/after spectrum plot window.")


def main():
    file_path = "data/example.imzML"

    with MSDataManagerImzML(filepath=file_path) as data_manager:
        data_manager.load_head_data()

        pre_first_signal, pre_first_noise, pre_first_snr = calculate_snr_details(data_manager.ms[0])
        pre_last_signal, pre_last_noise, pre_last_snr = calculate_snr_details(data_manager.ms[-1])

        # Use flat_generator + flat numba noise reduction for faster compute path.
        processed_data_manager = MSDataManagerImzML(MassSpectrumSet(), temp_dir="temp")
        processed_data_manager.copy_meta(data_manager)

        for mz_flat, intensity_flat, lengths, coordinates in data_manager.flat_generator(
            batch_size=512,
            include_mz=True,
            max_threads=16,
        ):
            nr_flat = FlatPreprocess.noise_reduction_flat(
                mz_data=mz_flat,
                intensity=intensity_flat,
                lengths=lengths,
                method="gaussian_numba",
            )

            processed_data_manager.swap_flat_data_out2disk(
                mz_flat=mz_flat,
                intensity_flat=nr_flat.intensity,
                lengths=lengths,
                coordinates=coordinates,
            )

        processed_data_manager.close_writer()
        processed_data_manager.load_head_data()

        logger.info(
            "Noise reduction finished with flat_generator (method=gaussian_numba). "
            f"Processed spectra: {len(processed_data_manager.ms)}"
        )

        post_first_signal, post_first_noise, post_first_snr = calculate_snr_details(processed_data_manager.ms[0])
        post_last_signal, post_last_noise, post_last_snr = calculate_snr_details(processed_data_manager.ms[-1])

        log_snr_details("Before processing - First spectrum", pre_first_signal, pre_first_noise, pre_first_snr)
        log_snr_details("Before processing - Last spectrum", pre_last_signal, pre_last_noise, pre_last_snr)
        log_snr_details("After processing - First spectrum", post_first_signal, post_first_noise, post_first_snr)
        log_snr_details("After processing - Last spectrum", post_last_signal, post_last_noise, post_last_snr)

        logger.info(
            "SNR change - First spectrum: %.4f -> %.4f (delta=%.4f)",
            pre_first_snr,
            post_first_snr,
            post_first_snr - pre_first_snr,
        )
        logger.info(
            "SNR change - Last spectrum: %.4f -> %.4f (delta=%.4f)",
            pre_last_snr,
            post_last_snr,
            post_last_snr - pre_last_snr,
        )

        plot_before_after(
            before_spectrum=data_manager.ms[0],
            after_spectrum=processed_data_manager.ms[0],
        )

        processed_data_manager.close()

if __name__ == "__main__":
    main()
