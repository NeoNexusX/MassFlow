"""Correctness tests for pixel-major MassFlow MSI Zarr stores."""

from __future__ import annotations

import numpy as np

from massflow.data_manager import MSDataManagerZarr
from massflow.msi_zarr import (
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
            reader.coordinates[:], [[1, 1, 1], [2, 1, 1]]
        )
        mz, intensity = reader.read_spectrum(1)
        np.testing.assert_array_equal(mz, [200.0, 201.0])
        np.testing.assert_array_equal(intensity, [2.0, 3.0])
    finally:
        reader.close()


def test_pixel_major_preserves_duplicate_coordinates_and_z(tmp_path):
    path = tmp_path / "duplicate_coordinates.zarr"
    writer = MSIZarrWriter(
        path,
        row_axis="pixel",
        encoding=ZARR_MODE_CONTINUOUS,
        metadata=None,
        width=2,
        height=2,
    )
    writer.write_flat_batch(
        intensity_flat=np.array([1.0, 10.0, 2.0, 20.0], dtype=np.float32),
        lengths=np.array([2, 2], dtype=np.uint32),
        coordinates=np.array([[1, 1, 7], [1, 1, 9]], dtype=np.uint32),
        mz_axis=np.array([100.0, 200.0]),
    )
    writer.finalize()
    writer.close()

    reader = MSIZarrReader(path)
    try:
        np.testing.assert_array_equal(
            reader.coordinates[:],
            [[1, 1, 7], [1, 1, 9]],
        )
        assert reader.n_rows == 2
        np.testing.assert_array_equal(reader.read_row(0), [1.0, 10.0])
        np.testing.assert_array_equal(reader.read_row(1), [2.0, 20.0])
        reader.validate_full()
    finally:
        reader.close()


def test_zarr_to_zarr_preserves_one_based_subset_coordinates(tmp_path):
    source_path = tmp_path / "source.zarr"
    writer = MSIZarrWriter(
        source_path,
        row_axis="pixel",
        encoding=ZARR_MODE_CONTINUOUS,
        metadata=None,
        width=3,
        height=2,
    )
    writer.write_flat_batch(
        intensity_flat=np.array([1.0, 10.0], dtype=np.float32),
        lengths=np.array([2], dtype=np.uint32),
        coordinates=np.array([[3, 2, 1]], dtype=np.uint32),
        mz_axis=np.array([100.0, 200.0]),
    )
    writer.finalize()
    writer.close()

    source = MSDataManagerZarr(filepath=source_path, mode="r")
    source.load_head_data()
    output = MSDataManagerZarr(
        filepath=tmp_path / "copy.zarr",
        mode="w",
        storage_mode=source.storage_mode,
    )
    output.copy_meta(source)
    for mz_axis, intensity, lengths, coordinates in source.flat_generator(
        batch_size=1,
        include_mz=True,
    ):
        output.write_flat_batch(intensity, lengths, coordinates, mz_axis)
    output.finalize()
    try:
        copied = MSIZarrReader(output.zarr_filepath)
        try:
            np.testing.assert_array_equal(copied.coordinates[:], [[3, 2, 1]])
        finally:
            copied.close()
    finally:
        output.close()
        source.close()
