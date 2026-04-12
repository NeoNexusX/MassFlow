from massflow.data_manager import MSDataManagerImzML
from massflow.tools import get_logger
from massflow.module import MassSpectrumSet
from massflow.preprocess.flat_pre_fun import FlatPreprocess
from massflow.preprocess.spectrum_pre_fun import SpectrumPreprocess
import numpy as np

logger = get_logger("massflow")


def calculate_snr_details(spectrum, method: str = "sd") -> tuple[float, float, float]:
    """Calculate signal level, noise level, and SNR for one spectrum."""
    intensity = spectrum.intensity
    signal_level = float(np.percentile(intensity, 95))
    noise = float(np.mean(SpectrumPreprocess.noise_estimation_spectrum(spectrum, method=method)))
    snr = float(signal_level / noise) if noise > 0 else float("inf")
    return signal_level, noise, snr


def log_snr_details(tag: str, signal_level: float, noise: float, snr: float) -> None:
    """Log signal level, noise level, and SNR for one spectrum."""
    logger.info(f"[{tag}] signal_level(95th)={signal_level:.4f}, noise={noise:.4f}, SNR={snr:.4f}")


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

        processed_data_manager.close()

if __name__ == "__main__":
    main()
