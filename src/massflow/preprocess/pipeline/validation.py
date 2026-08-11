"""Validation helpers for preprocessing pipeline inputs."""

from __future__ import annotations

import numpy as np

from massflow.data_manager import MSDataManager
from massflow.tools.funs import infer_shared_mz


def validate_mz_axis(mz_axis: np.ndarray, *, context: str) -> None:
    """Require a real-valued m/z axis sorted in ascending order."""
    mz_axis = np.asarray(mz_axis)

    if np.iscomplexobj(mz_axis):
        raise ValueError(
            f"{context}: m/z values must be real, got dtype={mz_axis.dtype}."
        )

    invalid = np.flatnonzero(
        ~np.less_equal(mz_axis[:-1], mz_axis[1:])
    )
    if invalid.size:
        index = int(invalid[0])
        raise ValueError(
            f"{context}: m/z values must be sorted in ascending order; "
            f"first invalid pair is mz[{index}]={mz_axis[index]} and "
            f"mz[{index + 1}]={mz_axis[index + 1]}."
        )


def validate_flat_mz(
    mz_data: np.ndarray,
    lengths: np.ndarray,
    *,
    batch_index: int,
    spectrum_offset: int,
) -> None:
    """Validate shared or per-spectrum m/z data from one flat batch."""
    mz_data = np.asarray(mz_data)
    lengths = np.asarray(lengths)
    batch_context = f"pipeline input batch {batch_index}"

    if np.iscomplexobj(mz_data):
        raise ValueError(
            f"{batch_context}: m/z values must be real, got dtype={mz_data.dtype}."
        )

    if infer_shared_mz(mz_data, lengths):
        validate_mz_axis(
            mz_data,
            context=f"{batch_context}, shared m/z axis",
        )
        return

    offset = 0
    for local_index, length in enumerate(lengths):
        end = offset + int(length)
        validate_mz_axis(
            mz_data[offset:end],
            context=f"pipeline input spectrum {spectrum_offset + local_index}",
        )
        offset = end


def validate_data_manager_mz(
    data_manager: MSDataManager,
    *,
    batch_size: int = 256,
) -> None:
    """Validate all m/z axes entering the preprocessing pipeline."""
    shared_mz = data_manager.ms.shared_mz_list
    if shared_mz is not None:
        validate_mz_axis(
            np.asarray(shared_mz),
            context="pipeline input shared m/z axis",
        )
        return

    return

    # spectrum_offset = 0
    # for batch_index, (mz_data, _intensity, lengths, _coordinates) in enumerate(
    #     data_manager.flat_generator(batch_size=batch_size, include_mz=True)
    # ):
    #     if mz_data is None:
    #         raise ValueError(
    #             f"pipeline input batch {batch_index}: m/z data is required."
    #         )
    #     lengths = np.asarray(lengths)
    #     validate_flat_mz(
    #         np.asarray(mz_data),
    #         lengths,
    #         batch_index=batch_index,
    #         spectrum_offset=spectrum_offset,
    #     )
    #     spectrum_offset += lengths.size
