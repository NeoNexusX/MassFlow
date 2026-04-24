from __future__ import annotations

import numpy as np

from massflow.tools.funs import prepare_flat_inputs, infer_shared_mz
from massflow.preprocess.helper.est_noise_helper import estimator
from massflow.preprocess.numba import peak_pick_numba_parallel as compute
from massflow.tools.logger import get_logger

logger = get_logger("massflow.peak_pick_parallel")


def peak_picker(
    mz_data: np.ndarray,
    intensity: np.ndarray,
    lengths: np.ndarray,
    width: int = 5,
    method: str = "quantile",
    snr: float = 2.0,
    return_type: str = "height",
    prominence: float | None = None,
    relheight: float | None = None,
    nbins: int = 1,
    overlap: float = 0.5,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Parallel peak picking for mass spectrometry data.
    """
    width = 1 if width < 1 else width

    mz_arr, intensity_arr, lengths_arr = prepare_flat_inputs(mz_data, intensity, lengths)

    is_shared_mz = infer_shared_mz(mz_arr, lengths_arr)

    # Execute everything in the numba parallel core: local maxima detection,
    # conditional filtering, and flat output generation.
    use_prominence = prominence is not None
    use_relheight = relheight is not None
    prominence = prominence if use_prominence else 0.0
    relheight = relheight if use_relheight else 0.0
    use_snr = np.isfinite(snr) and snr > 0.0

    noise_arr = np.empty(0, dtype=np.float64)

    if use_snr:
        noise_arr = estimator(
            intensity=intensity_arr,
            indexes=None,
            lengths=lengths_arr,
            nbins=nbins,
            overlap=overlap,
            dynamic=False,
            method=method,
            denoise_method="gaussian_numba",
        )

        if noise_arr.size != intensity_arr.size:
            raise ValueError("Noise estimator output size must match intensity_arr size.")

    return compute.peak_pick_flat_parallel_core_jit(
        mz_arr=mz_arr,
        intensity_arr=intensity_arr,
        lengths_arr=lengths_arr,
        width=width,
        use_prominence=use_prominence,
        prominence=prominence,
        use_relheight=use_relheight,
        relheight=relheight,
        use_snr=use_snr,
        snr=snr,
        return_area=(return_type == "area"),
        is_shared_mz=is_shared_mz,
        noise_arr=noise_arr,
    )
