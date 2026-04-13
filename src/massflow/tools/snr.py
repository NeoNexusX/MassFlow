from typing import Any

import numpy as np

from massflow.tools.logger import get_logger

logger = get_logger("massflow")


def calculate_snr_details(spectrum: Any, method: str = "sd") -> tuple[float, float, float]:
    """Calculate signal level, noise level, and SNR for one spectrum."""
    # Lazy import prevents circular import during package initialization.
    from massflow.preprocess.spectrum_pre_fun import SpectrumPreprocess

    intensity = spectrum.intensity
    signal_level = float(np.percentile(intensity, 95))
    noise = float(np.mean(SpectrumPreprocess.noise_estimation_spectrum(spectrum, method=method)))
    snr = float(signal_level / noise) if noise > 0 else float("inf")
    return signal_level, noise, snr


def log_snr_details(tag: str, signal_level: float, noise: float, snr: float) -> None:
    """Log signal level, noise level, and SNR for one spectrum."""
    logger.info(f"[{tag}] signal_level(95th)={signal_level:.4f}, noise={noise:.4f}, SNR={snr:.4f}")
