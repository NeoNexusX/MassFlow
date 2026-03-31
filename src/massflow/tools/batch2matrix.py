import numpy as np
from massflow.module import SpectrumImzML

def batch_to_matrix(batch):
    """convert a batch of spectra to a 2D mz matrix and intensity matrix."""

    lengths = np.array([len(s) for s in batch], dtype=np.int32)

    # 判定连续性 (Continuity Check)
    is_continuous = False
    if lengths[0] == lengths[-1] and np.all(lengths == lengths[0]):
        is_continuous = True

    if is_continuous:
        mz_matrix = np.vstack([s.mz_list for s in batch])
        intensity_matrix = np.vstack([s.intensity for s in batch])
    else:
        batch_size = len(batch)
        max_len = np.max(lengths)
        mz_matrix = np.zeros((batch_size, max_len), dtype=np.float64)
        intensity_matrix = np.zeros((batch_size, max_len), dtype=np.float32)
        for i, spec in enumerate(batch):
            if lengths[i] > 0:
                mz_matrix[i, :lengths[i]] = spec.mz_list
                intensity_matrix[i, :lengths[i]] = spec.intensity

    return mz_matrix, intensity_matrix, lengths

def matrix_to_batch(batch, mz_matrix, intensity_matrix, lengths):
    """convert a 2D mz matrix and intensity matrix back to a batch of spectra."""
    spectra: list[SpectrumImzML] = []
    shared_mz = None
    if mz_matrix.ndim == 1:
        shared_mz = mz_matrix

    for i, spectrum in enumerate(batch):
        valid_len = lengths[i]
        # Trim result back to original length
        mz = mz_matrix[i, :valid_len] if shared_mz is None else shared_mz
        intensity = intensity_matrix[i, :valid_len]

        spectra.append(
            SpectrumImzML(
                coordinates=spectrum.coordinate,
                mz_list=mz,
                intensity=intensity,
            )
        )
    return spectra
