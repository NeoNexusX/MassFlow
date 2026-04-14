from typing import Any

import numpy as np

from massflow.tools.logger import get_logger

logger = get_logger("massflow.tools.snr")


def calculate_snr_details(spectrum: Any, method: str = "sd") -> tuple[float, float, float]:
    """Calculate signal level, noise level, and SNR for one spectrum."""
    # Lazy import prevents circular import during package initialization.
    from massflow.preprocess.spectrum_pre_fun import SpectrumPreprocess

    intensity = np.asarray(spectrum.intensity)
    if intensity.size == 0 or not np.isfinite(intensity).any():
        msg = "Cannot calculate SNR details: spectrum intensity is empty or has no finite values."
        raise ValueError(msg)

    signal_level = float(np.percentile(intensity, 95))
    if not np.isfinite(signal_level):
        msg = "Cannot calculate SNR details: computed signal level is not finite."
        raise ValueError(msg)

    noise_values = np.asarray(SpectrumPreprocess.noise_estimation_spectrum(spectrum, method=method))
    if noise_values.size == 0 or not np.isfinite(noise_values).any():
        msg = "Cannot calculate SNR details: noise estimate is empty or has no finite values."
        raise ValueError(msg)

    noise = float(np.mean(noise_values))
    if not np.isfinite(noise):
        msg = "Cannot calculate SNR details: computed noise level is not finite."
        raise ValueError(msg)

    snr = float(signal_level / noise) if noise > 0 else float("inf")
    return signal_level, noise, snr


def log_snr_details(tag: str, signal_level: float, noise: float, snr: float) -> None:
    """Log signal level, noise level, and SNR for one spectrum."""
    logger.info(f"[{tag}] signal_level(95th)={signal_level:.4f}, noise={noise:.4f}, SNR={snr:.4f}")
