"""Correctness tests for pixel-major MassFlow MSI Zarr stores."""

from __future__ import annotations

import numpy as np

from massflow.data_manager import MSDataManagerZarr
from massflow.tools.zarr_io import (
    MSIZarrReader,
    MSIZarrWriter,
    ZARR_MODE_CONTINUOUS,
    ZARR_MODE_PROCESSED,
)


def test_continuous_pixel_major_roundtrip_and_add_flush(tmp_path):
    path = tmp_path / "continuous.zarr"
    writer = MSIZarrWriter(
        path,
        row_axis="pixel",
        encoding=ZARR_MODE_CONTINUOUS,
        metadata=None,
        width=3,
        height=2,
        coordinate_base=0,
    )
    writer.add_spectrum(
        np.array([1.0, 10.0], dtype=np.float32),
        np.array([2, 1, 0], dtype=np.uint32),
        np.array([100.0, 200.0]),
    )
    writer.finalize()
    writer.close()

    reader = MSIZarrReader(path)
    try:
        assert reader.row_axis == "pixel"
        assert reader.encoding == ZARR_MODE_CONTINUOUS
        np.testing.assert_array_equal(reader.coordinates[:], [[2, 1, 0]])
        mz, intensity = reader.read_spectrum(0)
        np.testing.assert_array_equal(mz, [100.0, 200.0])
        np.testing.assert_array_equal(intensity, [1.0, 10.0])
    finally:
        reader.close()


def test_processed_pixel_major_roundtrip(tmp_path):
    path = tmp_path / "processed.zarr"
    writer = MSIZarrWriter(
        path,
        row_axis="pixel",
        encoding=ZARR_MODE_PROCESSED,
        metadata=None,
        width=2,
        height=1,
        coordinate_base=1,
    )
    writer.write_flat_batch(
        intensity_flat=np.array([1.0, 2.0, 3.0], dtype=np.float32),
        lengths=np.array([1, 2], dtype=np.uint32),
        coordinates=np.array([[1, 1, 1], [2, 1, 1]], dtype=np.uint32),
        mz_axis=np.array([100.0, 200.0, 201.0]),
    )
    writer.finalize()
    writer.close()

    reader = MSIZarrReader(path)
    try:
        assert reader.encoding == ZARR_MODE_PROCESSED
        np.testing.assert_array_equal(
            reader.coordinates[:], [[0, 0, 0], [1, 0, 0]]
        )
        mz, intensity = reader.read_spectrum(1)
        np.testing.assert_array_equal(mz, [200.0, 201.0])
        np.testing.assert_array_equal(intensity, [2.0, 3.0])
    finally:
        reader.close()


def test_zarr_to_zarr_preserves_zero_based_subset_coordinates(tmp_path):
    source_path = tmp_path / "source.zarr"
    writer = MSIZarrWriter(
        source_path,
        row_axis="pixel",
        encoding=ZARR_MODE_CONTINUOUS,
        metadata=None,
        width=3,
        height=2,
        coordinate_base=0,
    )
    writer.write_flat_batch(
        intensity_flat=np.array([1.0, 10.0], dtype=np.float32),
        lengths=np.array([2], dtype=np.uint32),
        coordinates=np.array([[2, 1, 0]], dtype=np.uint32),
        mz_axis=np.array([100.0, 200.0]),
    )
    writer.finalize()
    writer.close()

    source = MSDataManagerZarr(filepath=source_path, mode="r")
    source.load_head_data()
    output = MSDataManagerZarr.from_data_manager(
        source,
        tmp_path / "copy.zarr",
        batch_size=1,
    )
    try:
        copied = MSIZarrReader(output.zarr_filepath)
        try:
            np.testing.assert_array_equal(copied.coordinates[:], [[2, 1, 0]])
        finally:
            copied.close()
    finally:
        output.close()
        source.close()
