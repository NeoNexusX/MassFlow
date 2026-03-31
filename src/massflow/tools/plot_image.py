"""Image plotting helpers for MSI data."""

from __future__ import annotations

from typing import Optional

import matplotlib.pyplot as plt
import numpy as np

from massflow.tools.logger import get_logger


logger = get_logger("massflow.tools.plot_image")


def plot_single_ion_image(
    processed_data_manager,
    target_mz: float,
    mz_tolerance: float = 0.05,
    batch_size: int = 256,
    figsize: tuple[float, float] = (6, 5),
    title: Optional[str] = None,
):
    """Plot a single-ion image from a processed data manager.

    The color upper bound is fixed at the 95th percentile, the colormap is
    inferno, and a colorbar is always attached.

    Args:
        processed_data_manager: Data manager containing processed MSI spectra.
        target_mz: Target m/z value for the ion image.
        mz_tolerance: Absolute m/z tolerance around ``target_mz``.
        batch_size: Batch size used by ``matrix_generator``.
        figsize: Matplotlib figure size.
        title: Optional title. If omitted, an automatic title is used.

    Returns:
        tuple[np.ndarray, matplotlib.figure.Figure, matplotlib.axes.Axes]:
            Image matrix, figure, and axes.
    """
    meta = getattr(processed_data_manager.ms, "meta", None)
    if meta is None:
        raise ValueError("processed_data_manager.ms.meta is required for image plotting.")

    width = int(meta.max_count_of_pixels_x)
    height = int(meta.max_count_of_pixels_y)
    if width <= 0 or height <= 0:
        raise ValueError("Invalid image size from metadata.")

    img_data = np.zeros((height, width), dtype=np.float32)

    for mz_data, intensity_matrix, lengths, coords in processed_data_manager.matrix_generator(
        batch_size=batch_size,
        include_mz=True,
    ):
        if intensity_matrix.size == 0:
            continue
        if mz_data is None:
            raise ValueError("matrix_generator returned no m/z data; include_mz=True is required.")

        n = intensity_matrix.shape[0]
        ion_vals = np.zeros(n, dtype=np.float64)

        # Continuous mode uses one shared m/z axis for all spectra in this batch.
        if np.ndim(mz_data) == 1:
            mz_mask = np.abs(np.asarray(mz_data) - target_mz) <= mz_tolerance
            if np.any(mz_mask):
                ion_vals = intensity_matrix[:, mz_mask].sum(axis=1, dtype=np.float64)
        else:
            mz_matrix = np.asarray(mz_data)
            for i in range(n):
                valid_len = int(lengths[i])
                if valid_len <= 0:
                    continue
                mz_row = mz_matrix[i, :valid_len]
                mz_mask = np.abs(mz_row - target_mz) <= mz_tolerance
                if np.any(mz_mask):
                    ion_vals[i] = intensity_matrix[i, :valid_len][mz_mask].sum(dtype=np.float64)

        x_idx = coords[:, 0].astype(np.int64)
        y_idx = coords[:, 1].astype(np.int64)

        # Handle common imzML 1-based coordinates.
        if x_idx.min() >= 1 and y_idx.min() >= 1:
            x_idx = x_idx - 1
            y_idx = y_idx - 1

        valid = (x_idx >= 0) & (x_idx < width) & (y_idx >= 0) & (y_idx < height)
        img_data[y_idx[valid], x_idx[valid]] += ion_vals[valid]

    non_zero = img_data[img_data > 0]
    vmax_val = np.percentile(non_zero, 95) if non_zero.size > 0 else np.percentile(img_data, 95)

    fig, ax = plt.subplots(figsize=figsize)
    im = ax.imshow(img_data, cmap="inferno", vmax=vmax_val)
    ax.set_title(title or f"Single Ion Image (m/z={target_mz:.4f}, tol={mz_tolerance:g})")
    ax.set_xlabel("X")
    ax.set_ylabel("Y")

    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("Intensity")

    plt.tight_layout()
    plt.show()

    logger.info(
        "Single ion image plotted: target_mz=%s, tolerance=%s, shape=(%s, %s)",
        target_mz,
        mz_tolerance,
        height,
        width,
    )
    return img_data, fig, ax


def plot_tic_image(
    processed_data_manager,
    batch_size: int = 256,
    figsize: tuple[float, float] = (6, 5),
    title: Optional[str] = None,
):
    """Plot a total ion current (TIC) image from a processed data manager.

    The color upper bound is fixed at the 95th percentile, the colormap is
    inferno, and a colorbar is always attached.

    Args:
        processed_data_manager: Data manager containing processed MSI spectra.
        batch_size: Batch size used by ``matrix_generator``.
        figsize: Matplotlib figure size.
        title: Optional title. If omitted, an automatic title is used.

    Returns:
        tuple[np.ndarray, matplotlib.figure.Figure, matplotlib.axes.Axes]:
            Image matrix, figure, and axes.
    """
    meta = getattr(processed_data_manager.ms, "meta", None)
    if meta is None:
        raise ValueError("processed_data_manager.ms.meta is required for image plotting.")

    width = int(meta.max_count_of_pixels_x)
    height = int(meta.max_count_of_pixels_y)
    if width <= 0 or height <= 0:
        raise ValueError("Invalid image size from metadata.")

    img_data = np.zeros((height, width), dtype=np.float32)

    for _, intensity_matrix, lengths, coords in processed_data_manager.matrix_generator(
        batch_size=batch_size,
        include_mz=False,
    ):
        if intensity_matrix.size == 0:
            continue

        n = intensity_matrix.shape[0]
        tic_vals = np.zeros(n, dtype=np.float64)

        if lengths is None:
            tic_vals = intensity_matrix.sum(axis=1, dtype=np.float64)
        else:
            for i in range(n):
                valid_len = int(lengths[i])
                if valid_len <= 0:
                    continue
                tic_vals[i] = intensity_matrix[i, :valid_len].sum(dtype=np.float64)

        x_idx = coords[:, 0].astype(np.int64)
        y_idx = coords[:, 1].astype(np.int64)

        # Handle common imzML 1-based coordinates.
        if x_idx.min() >= 1 and y_idx.min() >= 1:
            x_idx = x_idx - 1
            y_idx = y_idx - 1

        valid = (x_idx >= 0) & (x_idx < width) & (y_idx >= 0) & (y_idx < height)
        img_data[y_idx[valid], x_idx[valid]] += tic_vals[valid]

    non_zero = img_data[img_data > 0]
    vmax_val = np.percentile(non_zero, 95) if non_zero.size > 0 else np.percentile(img_data, 95)

    fig, ax = plt.subplots(figsize=figsize)
    im = ax.imshow(img_data, cmap="inferno", vmax=vmax_val)
    ax.set_title(title or "TIC Image")
    ax.set_xlabel("X")
    ax.set_ylabel("Y")

    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("Intensity")

    plt.tight_layout()
    plt.show()

    logger.info(
        "TIC image plotted: shape=(%s, %s)",
        height,
        width,
    )
    return img_data, fig, ax
