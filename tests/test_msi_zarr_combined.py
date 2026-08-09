"""Combined ion_image + spectra v1.1 MSI Zarr schema tests."""

# pylint: disable=redefined-outer-name,protected-access,too-many-lines

from __future__ import annotations

import shutil
from pathlib import Path
from typing import cast

import numpy as np
import pytest
import zarr

from massflow.data_manager import IonImageDataManagerZarr, MSDataManagerZarr
from massflow.msi_zarr import (
    MSIZarrReader,
    MSIZarrWriter,
    ZARR_MODE_CONTINUOUS,
    ZARR_MODE_PROCESSED,
)
from massflow.segmentation.region_mask import apply_region_mask

MZ_AXIS = np.array([100.0, 200.0])


def _write_combined(
    path: Path,
    tmp_path: Path,
    spectra_rows,
    *,
    mz: np.ndarray = MZ_AXIS,
    coords=None,
    width: int | None = None,
    height: int | None = None,
) -> IonImageDataManagerZarr:
    """Write one batch of pixel-major rows as a combined v1.1 store."""
    rows = np.asarray(spectra_rows, dtype=np.float32)
    if rows.ndim == 1:
        rows = rows.reshape(1, -1)
    n_pixels = rows.shape[0]
    if coords is None:
        coords = np.array(
            [[i + 1, 1, 1] for i in range(n_pixels)], dtype=np.int64
        )
    dm = IonImageDataManagerZarr(
        str(path),
        mode="w",
        width=width,
        height=height,
        temp_dir=str(tmp_path),
        intensity_dtype=np.float32,
        ion_block=16,
        include_spectra=True,
    )
    try:
        dm.write_flat_batch(
            intensity_flat=rows.ravel(),
            lengths=np.full(n_pixels, rows.shape[1], dtype=np.uint32),
            coordinates=coords,
            mz_axis=np.asarray(mz, dtype=np.float64),
        )
        dm.finalize()
        return dm
    except BaseException:
        dm.close()
        raise


def _copy_store(src: Path, dst: Path) -> None:
    shutil.copytree(src, dst)


def test_ion_only_v1_1_no_spectra(tmp_path):
    """include_spectra=False writes /ion_image only and keeps ion API values."""
    path = tmp_path / "ion_only.zarr"
    dm = IonImageDataManagerZarr(
        str(path),
        mode="w",
        width=2,
        height=1,
        temp_dir=str(tmp_path),
        intensity_dtype=np.float32,
        ion_block=16,
    )
    try:
        dm.write_flat_batch(
            intensity_flat=np.array([1, 2, 10, 20], dtype=np.float32),
            lengths=np.array([2, 2], dtype=np.uint32),
            coordinates=np.array([[1, 1, 1], [2, 1, 1]], dtype=np.uint32),
            mz_axis=np.array([100.0, 200.0]),
        )
        dm.finalize()
    finally:
        dm.close()

    root = zarr.open_group(str(path), mode="r")
    assert root.attrs["format_version"] == "1.1"
    assert "ion_image" in root
    assert "spectra" not in root
    assert "data" not in root
    reader = MSIZarrReader(path)
    try:
        assert reader.row_axis == "ion"
        np.testing.assert_array_equal(reader.read_ion_vector(0), [1.0, 2.0])
        np.testing.assert_array_equal(reader.read_ion_vector(1), [10.0, 20.0])
    finally:
        reader.close()


@pytest.mark.parametrize("encoding", [ZARR_MODE_CONTINUOUS, ZARR_MODE_PROCESSED])
def test_standalone_pixel_major_uses_spectra(tmp_path, encoding):
    """Standalone pixel-major writer emits v1.1 /spectra without /data."""
    path = tmp_path / f"spec_{encoding}.zarr"
    writer = MSIZarrWriter(
        path,
        row_axis="pixel",
        encoding=encoding,
        metadata=None,
        width=2,
        height=1,
    )
    if encoding == ZARR_MODE_CONTINUOUS:
        mz = np.array([100.0, 200.0])
    else:
        mz = np.array([100.0, 100.1, 200.0, 200.1])
    writer.write_flat_batch(
        intensity_flat=np.array([1.0, 10.0, 2.0, 20.0], dtype=np.float32),
        lengths=np.array([2, 2], dtype=np.uint32),
        coordinates=np.array([[1, 1, 1], [2, 1, 1]], dtype=np.uint32),
        mz_axis=mz,
    )
    writer.finalize()
    writer.close()

    root = zarr.open_group(str(path), mode="r")
    assert root.attrs["format_version"] == "1.1"
    assert "spectra" in root
    assert "ion_image" not in root
    assert "data" not in root
    reader = MSIZarrReader(path)
    try:
        assert reader.row_axis == "pixel"
        mz_read, intensity = reader.read_spectrum(1)
        np.testing.assert_array_equal(intensity, [2.0, 20.0])
        if encoding == ZARR_MODE_CONTINUOUS:
            np.testing.assert_array_equal(mz_read, [100.0, 200.0])
        else:
            np.testing.assert_array_equal(mz_read, [200.0, 200.1])
    finally:
        reader.close()


def test_combined_layout_and_known_matrix(tmp_path):
    """Both layouts coexist, share axes, and match the known 2x2 example."""
    path = tmp_path / "known.zarr"
    dm = _write_combined(path, tmp_path, [[1.0, 10.0], [2.0, 20.0]])
    dm.close()

    root = zarr.open_group(str(path), mode="r")
    assert dict(root.attrs)["format_version"] == "1.1"
    assert "ion_image" in root and "spectra" in root
    assert "data" not in root
    assert dict(root["ion_image"].attrs) == {
        "row_axis": "ion",
        "encoding": "continuous",
    }
    assert dict(root["spectra"].attrs) == {
        "row_axis": "pixel",
        "encoding": "continuous",
    }
    np.testing.assert_array_equal(
        cast(zarr.Array, root["spectra/intensity"])[:],
        [1.0, 10.0, 2.0, 20.0],
    )
    np.testing.assert_array_equal(
        cast(zarr.Array, root["spectra/offsets"])[:], [0, 2, 4]
    )
    np.testing.assert_array_equal(
        cast(zarr.Array, root["ion_image/intensity"])[:],
        [1.0, 2.0, 10.0, 20.0],
    )
    np.testing.assert_array_equal(
        cast(zarr.Array, root["ion_image/offsets"])[:], [0, 2, 4]
    )
    np.testing.assert_array_equal(
        cast(zarr.Array, root["axes/mz"])[:], [100.0, 200.0]
    )
    np.testing.assert_array_equal(
        cast(zarr.Array, root["axes/coordinates"])[:],
        [[1, 1, 1], [2, 1, 1]],
    )


def test_combined_multibatch_pixel_order_alignment(tmp_path):
    """Multiple batches keep /spectra and /ion_image pixel orders identical."""
    path = tmp_path / "multi.zarr"
    dm = IonImageDataManagerZarr(
        str(path),
        mode="w",
        width=3,
        height=1,
        temp_dir=str(tmp_path),
        intensity_dtype=np.float32,
        ion_block=16,
        include_spectra=True,
    )
    try:
        for coord, values in [
            ([1, 1, 1], [1.0, 10.0]),
            ([2, 1, 1], [3.0, 30.0]),
            ([3, 1, 1], [4.0, 40.0]),
        ]:
            dm.write_flat_batch(
                intensity_flat=np.asarray(values, dtype=np.float32),
                lengths=np.array([2], dtype=np.uint32),
                coordinates=np.asarray([coord], dtype=np.uint32),
                mz_axis=np.array([100.0, 200.0]),
            )
        dm.finalize()
    finally:
        dm.close()

    reader = MSIZarrReader(path)
    try:
        np.testing.assert_array_equal(reader.coordinates[:], [[1, 1, 1], [2, 1, 1], [3, 1, 1]])
        np.testing.assert_array_equal(reader.read_ion_vector(0), [1.0, 3.0, 4.0])
        np.testing.assert_array_equal(reader.read_ion_vector(1), [10.0, 30.0, 40.0])
        assert reader._spectra_group is not None
        np.testing.assert_array_equal(
            np.asarray(reader._spectra_group["intensity"][:]),
            [1.0, 10.0, 3.0, 30.0, 4.0, 40.0],
        )
        np.testing.assert_array_equal(
            np.asarray(reader._spectra_group["offsets"][:]), [0, 2, 4, 6]
        )
        reader.validate_full()
    finally:
        reader.close()


def test_combined_duplicate_coordinates_compact(tmp_path):
    """Duplicate coordinates compact both layouts and keep the last TIC."""
    path = tmp_path / "dup.zarr"
    dm = IonImageDataManagerZarr(
        str(path),
        mode="w",
        width=2,
        height=1,
        temp_dir=str(tmp_path),
        intensity_dtype=np.float32,
        ion_block=16,
        include_spectra=True,
    )
    try:
        dm.write_flat_batch(
            intensity_flat=np.array([1.0, 10.0], dtype=np.float32),
            lengths=np.array([2], dtype=np.uint32),
            coordinates=np.array([[1, 1, 1]], dtype=np.uint32),
            mz_axis=np.array([100.0, 200.0]),
        )
        dm.write_flat_batch(
            intensity_flat=np.array([2.0, 20.0], dtype=np.float32),
            lengths=np.array([2], dtype=np.uint32),
            coordinates=np.array([[1, 1, 1]], dtype=np.uint32),
            mz_axis=np.array([100.0, 200.0]),
        )
        dm.finalize()
    finally:
        dm.close()

    reader = MSIZarrReader(path)
    try:
        assert reader.n_pixels == 1
        np.testing.assert_array_equal(reader.coordinates[:], [[1, 1, 1]])
        np.testing.assert_array_equal(reader.read_ion_vector(0), [2.0])
        np.testing.assert_array_equal(reader.read_ion_vector(1), [20.0])
        assert reader._spectra_group is not None
        np.testing.assert_array_equal(
            np.asarray(reader._spectra_group["intensity"][:]),
            [2.0, 20.0],
        )
        np.testing.assert_array_equal(
            np.asarray(reader._spectra_group["offsets"][:]), [0, 2]
        )
        tic = reader.read_tic()
        assert tic is not None
        np.testing.assert_allclose(tic, [22.0])
        mean = reader.read_mean_spectrum()
        assert mean is not None
        np.testing.assert_allclose(mean, [2.0, 20.0])
        reader.validate_full()
    finally:
        reader.close()


def test_combined_holes_and_shape_inference(tmp_path):
    """Spatial holes and one-based width/height inference are preserved."""
    path = tmp_path / "holes.zarr"
    dm = _write_combined(
        path,
        tmp_path,
        [[1.0, 10.0], [2.0, 20.0]],
        coords=np.array([[1, 1, 1], [3, 1, 1]], dtype=np.int64),
        width=None,
        height=None,
    )
    try:
        assert dm.width == 3
        assert dm.height == 1
        np.testing.assert_array_equal(
            dm.get_ion_image(0), np.array([[1, 0, 2]], dtype=np.float32)
        )
    finally:
        dm.close()


def test_combined_reads_via_ion_major_api(tmp_path):
    """Reader and IonImageDataManagerZarr expose combined files as ion-major."""
    path = tmp_path / "combined.zarr"
    dm = _write_combined(path, tmp_path, [[1.0, 10.0], [2.0, 20.0]])
    dm.close()

    reader = MSIZarrReader(path)
    try:
        assert reader.row_axis == "ion"
        assert reader.encoding == ZARR_MODE_CONTINUOUS
        np.testing.assert_array_equal(reader.read_ion_vector(0), [1.0, 2.0])
        np.testing.assert_array_equal(
            reader.read_ion_image(0), np.array([[1, 2]], dtype=np.float32)
        )
        reader.validate_full()
    finally:
        reader.close()

    with IonImageDataManagerZarr(str(path), mode="r") as reopened:
        np.testing.assert_array_equal(reopened.get_ion_vector(0), [1.0, 2.0])
        np.testing.assert_array_equal(
            reopened.get_ion_image(0), np.array([[1, 2]], dtype=np.float32)
        )


def test_reader_rejects_v1_0_and_data_group(tmp_path):
    """Reader refuses format_version 1.0 and any store with a /data group."""
    path = tmp_path / "v1.zarr"
    root = zarr.open_group(str(path), mode="w")
    root.attrs.update(
        {
            "format": "massflow.msi_zarr",
            "format_version": "1.0",
            "write_state": "complete",
            "coordinate_order": ["x", "y", "z"],
            "coordinate_base": 1,
            "spatial_shape": [1, 2],
        }
    )
    root.create_group("metadata")
    axes = root.create_group("axes")
    axes.create_array(
        "coordinates", data=np.array([[1, 1, 1], [2, 1, 1]], dtype="<u4")
    )
    mz = axes.create_array("mz", data=np.array([100.0, 200.0]))
    mz.attrs["unit"] = "m/z"
    data = root.create_group("data")
    data.attrs.update({"row_axis": "ion", "encoding": "continuous"})
    data.create_array(
        "intensity", data=np.array([1.0, 2.0, 10.0, 20.0], dtype=np.float32)
    )
    data.create_array("offsets", data=np.array([0, 2, 4], dtype="<i8"))
    root["metadata"].attrs["continuous"] = True
    with pytest.raises(ValueError, match="Unsupported.*version"):
        MSIZarrReader(path)

    path2 = tmp_path / "data_group.zarr"
    root2 = zarr.open_group(str(path2), mode="w")
    root2.attrs.update(
        {
            "format": "massflow.msi_zarr",
            "format_version": "1.1",
            "write_state": "complete",
            "coordinate_order": ["x", "y", "z"],
            "coordinate_base": 1,
            "spatial_shape": [1, 2],
        }
    )
    root2.create_group("metadata")
    root2.create_group("axes")
    root2.create_group("data")
    with pytest.raises(ValueError, match="data/ group"):
        MSIZarrReader(path2)


def test_combined_corruption_rejected(tmp_path):
    """Missing/incorrect /spectra contents reject the store on open."""
    base = tmp_path / "base.zarr"
    dm = _write_combined(base, tmp_path, [[1.0, 10.0], [2.0, 20.0]])
    dm.close()

    # Missing spectra/intensity
    path = tmp_path / "missing_intensity.zarr"
    _copy_store(base, path)
    root = zarr.open_group(str(path), mode="a")
    del root["spectra/intensity"]
    with pytest.raises(ValueError, match="spectra/intensity and spectra/offsets"):
        MSIZarrReader(path)

    # Wrong spectra row_axis
    path = tmp_path / "wrong_attrs.zarr"
    _copy_store(base, path)
    root = zarr.open_group(str(path), mode="a")
    root["spectra"].attrs["row_axis"] = "ion"
    with pytest.raises(ValueError, match="spectra attrs.row_axis"):
        MSIZarrReader(path)

    # Broken spectra offsets boundary
    path = tmp_path / "bad_offsets.zarr"
    _copy_store(base, path)
    root = zarr.open_group(str(path), mode="a")
    cast(zarr.Array, root["spectra/offsets"])[-1] = 99
    with pytest.raises(ValueError, match="spectra/offsets boundaries"):
        MSIZarrReader(path)

    # spectra row count mismatches axes/coordinates
    path = tmp_path / "pixel_mismatch.zarr"
    _copy_store(base, path)
    root = zarr.open_group(str(path), mode="a")
    offsets = cast(zarr.Array, root["spectra/offsets"])
    offsets.resize((offsets.shape[0] + 1,))
    offsets[-1] = int(offsets[-2])
    with pytest.raises(ValueError, match="spectra row count"):
        MSIZarrReader(path)


def test_include_spectra_rejected_for_pixel_major(tmp_path):
    """include_spectra=True is only valid for continuous ion-major output."""
    with pytest.raises(ValueError, match="include_spectra"):
        MSIZarrWriter(
            tmp_path / "pixel_continuous.zarr",
            row_axis="pixel",
            encoding=ZARR_MODE_CONTINUOUS,
            metadata=None,
            width=2,
            height=1,
            include_spectra=True,
        )
    with pytest.raises(ValueError, match="include_spectra"):
        MSIZarrWriter(
            tmp_path / "pixel_processed.zarr",
            row_axis="pixel",
            encoding=ZARR_MODE_PROCESSED,
            metadata=None,
            width=2,
            height=1,
            include_spectra=True,
        )
    with pytest.raises(ValueError, match="include_spectra"):
        MSIZarrWriter(
            tmp_path / "ion_auto.zarr",
            row_axis="ion",
            encoding=ZARR_MODE_PROCESSED,
            metadata=None,
            width=2,
            height=1,
            include_spectra=True,
        )


def test_combined_stats_generated_by_default(tmp_path):
    """include_spectra=True produces TIC and mean spectrum without extra flags."""
    path = tmp_path / "stats.zarr"
    dm = _write_combined(path, tmp_path, [[1.0, 10.0], [2.0, 20.0]])
    try:
        mean = dm.mean_spectrum
        assert mean is not None
        np.testing.assert_allclose(mean, [1.5, 15.0])
    finally:
        dm.close()

    reader = MSIZarrReader(path)
    try:
        tic = reader.read_tic()
        assert tic is not None
        np.testing.assert_allclose(tic, [11.0, 22.0])
        mean = reader.read_mean_spectrum()
        assert mean is not None
        np.testing.assert_allclose(mean, [1.5, 15.0])
        reader.validate_full()
    finally:
        reader.close()


def test_combined_stats_attrs_shared_group(tmp_path):
    """Both stats arrays carry exact dtype, shape and attrs in one stats group."""
    path = tmp_path / "attrs.zarr"
    dm = _write_combined(path, tmp_path, [[1.0, 10.0], [2.0, 20.0]])
    dm.close()

    root = zarr.open_group(str(path), mode="r")
    stats = cast(zarr.Group, root["stats"])
    tic = cast(zarr.Array, stats["tic"])
    mean = cast(zarr.Array, stats["mean_spectrum"])
    assert tic.shape == (2,)
    assert tic.dtype == np.float64
    assert tic.attrs["source"] == "/spectra"
    assert tic.attrs["reduced_axis"] == "mz"
    assert mean.shape == (2,)
    assert mean.dtype == np.float64
    assert mean.attrs["source"] == "/ion_image"
    assert mean.attrs["reduced_axis"] == "pixel"
    assert mean.attrs["denominator"] == "sampled_pixels"


def test_single_layout_stats_optional(tmp_path):
    """Ion-only and standalone spectra-only keep stats switches optional."""
    ion_path = tmp_path / "ion.zarr"
    dm = IonImageDataManagerZarr(
        str(ion_path),
        mode="w",
        width=2,
        height=1,
        temp_dir=str(tmp_path),
        intensity_dtype=np.float32,
        ion_block=16,
    )
    try:
        dm.write_flat_batch(
            intensity_flat=np.array([1.0, 2.0], dtype=np.float32),
            lengths=np.array([2], dtype=np.uint32),
            coordinates=np.array([[1, 1, 1]], dtype=np.uint32),
            mz_axis=np.array([100.0, 200.0]),
        )
        dm.finalize()
    finally:
        dm.close()
    root = zarr.open_group(str(ion_path), mode="r")
    assert "stats" not in root
    assert dm.mean_spectrum is None

    spec_path = tmp_path / "spec.zarr"
    writer = MSIZarrWriter(
        spec_path,
        row_axis="pixel",
        encoding=ZARR_MODE_CONTINUOUS,
        metadata=None,
        width=1,
        height=1,
    )
    writer.write_flat_batch(
        intensity_flat=np.array([1.0, 2.0], dtype=np.float32),
        lengths=np.array([2], dtype=np.uint32),
        coordinates=np.array([[1, 1, 1]], dtype=np.uint32),
        mz_axis=np.array([100.0, 200.0]),
    )
    writer.finalize()
    writer.close()
    root = zarr.open_group(str(spec_path), mode="r")
    assert "stats" not in root


def test_region_mask_writes_v1_1_ion_image(tmp_path):
    """region_mask output uses v1.1 /ion_image and reopens with the new reader."""
    source = IonImageDataManagerZarr(
        str(tmp_path / "source.zarr"),
        mode="w",
        width=2,
        height=1,
        temp_dir=str(tmp_path),
        intensity_dtype=np.float32,
        ion_block=16,
    )
    try:
        source.write_flat_batch(
            intensity_flat=np.array([1.0, 2.0, 10.0, 20.0], dtype=np.float32),
            lengths=np.array([2, 2], dtype=np.uint32),
            coordinates=np.array([[1, 1, 1], [2, 1, 1]], dtype=np.uint32),
            mz_axis=np.array([100.0, 200.0]),
        )
        source.finalize()
    finally:
        source.close()

    output: IonImageDataManagerZarr | None = None
    try:
        mask = np.array([[False, True]], dtype=bool)
        output = cast(
            IonImageDataManagerZarr,
            apply_region_mask(source, mask, tmp_path / "masked.zarr"),
        )
        reader = MSIZarrReader(output.filepath)
        try:
            assert reader.root.attrs["format_version"] == "1.1"
            assert "ion_image" in reader.root
            assert "spectra" not in reader.root
            assert "data" not in reader.root
            np.testing.assert_array_equal(reader.coordinates[:], [[2, 1, 1]])
            np.testing.assert_array_equal(
                reader.read_ion_vector(0), np.array([2.0], dtype=np.float32)
            )
            reader.validate_full()
        finally:
            reader.close()
    finally:
        if output is not None:
            output.close()


def test_combined_add_spectrum_buffer_flush(tmp_path):
    """add_spectrum() flushed through the buffer produces a combined store."""
    path = tmp_path / "buffer.zarr"
    writer = MSIZarrWriter(
        path,
        row_axis="ion",
        encoding=ZARR_MODE_CONTINUOUS,
        metadata=None,
        width=2,
        height=1,
        include_spectra=True,
        add_buffer_size=2,
    )
    writer.add_spectrum(
        np.array([1.0, 10.0], dtype=np.float32),
        np.array([1, 1, 1], dtype=np.uint32),
        np.array([100.0, 200.0]),
    )
    writer.add_spectrum(
        np.array([2.0, 20.0], dtype=np.float32),
        np.array([2, 1, 1], dtype=np.uint32),
        np.array([100.0, 200.0]),
    )
    writer.add_spectrum(
        np.array([3.0, 30.0], dtype=np.float32),
        np.array([2, 1, 1], dtype=np.uint32),
        np.array([100.0, 200.0]),
    )
    writer.finalize()
    writer.close()

    reader = MSIZarrReader(path)
    try:
        assert reader.n_pixels == 2  # duplicate (2,1,1) keeps the last spectrum
        assert reader._spectra_group is not None
        np.testing.assert_array_equal(
            np.asarray(reader._spectra_group["intensity"][:]),
            [1.0, 10.0, 3.0, 30.0],
        )
        np.testing.assert_array_equal(
            np.asarray(reader._spectra_group["offsets"][:]), [0, 2, 4]
        )
        np.testing.assert_array_equal(reader.read_ion_vector(1), [10.0, 30.0])
        tic = reader.read_tic()
        assert tic is not None
        np.testing.assert_allclose(tic, [11.0, 33.0])
        reader.validate_full()
    finally:
        reader.close()


def test_factories_forward_include_spectra(tmp_path):
    """Both IonImageDataManagerZarr factories forward include_spectra=True."""
    src_path = tmp_path / "src.zarr"
    writer = MSIZarrWriter(
        src_path,
        row_axis="pixel",
        encoding=ZARR_MODE_CONTINUOUS,
        metadata=None,
        width=2,
        height=1,
    )
    writer.write_flat_batch(
        intensity_flat=np.array([1.0, 10.0, 2.0, 20.0], dtype=np.float32),
        lengths=np.array([2, 2], dtype=np.uint32),
        coordinates=np.array([[1, 1, 1], [2, 1, 1]], dtype=np.uint32),
        mz_axis=np.array([100.0, 200.0]),
    )
    writer.finalize()
    writer.close()

    source = MSDataManagerZarr(filepath=src_path, mode="r")
    source.load_head_data()
    outputs: list[Path] = []
    try:
        dm = IonImageDataManagerZarr.from_data_manager(
            source, tmp_path / "from_dm.zarr", include_spectra=True
        )
        try:
            for mz, intensity, lengths, coords in source.flat_generator(
                batch_size=1, include_mz=True
            ):
                dm.write_flat_batch(intensity, lengths, coords, mz)
            dm.finalize()
            outputs.append(dm.filepath)
        finally:
            dm.close()

        dm2 = IonImageDataManagerZarr.from_flat_output_contract(
            source,
            tmp_path / "from_contract.zarr",
            shared_mz_axis=np.array([100.0, 200.0]),
            include_spectra=True,
        )
        try:
            for mz, intensity, lengths, coords in source.flat_generator(
                batch_size=1, include_mz=True
            ):
                dm2.write_flat_batch(intensity, lengths, coords, mz)
            dm2.finalize()
            outputs.append(dm2.filepath)
        finally:
            dm2.close()
    finally:
        source.close()

    for out in outputs:
        root = zarr.open_group(str(out), mode="r")
        assert root.attrs["format_version"] == "1.1"
        assert "ion_image" in root and "spectra" in root
        assert "data" not in root
        assert "stats/tic" in root and "stats/mean_spectrum" in root


def test_reader_rejects_dangling_stats_sources(tmp_path):
    """Deleting a layout group leaves a dangling stats source the reader rejects."""
    base = tmp_path / "base.zarr"
    dm = _write_combined(base, tmp_path, [[1.0, 10.0], [2.0, 20.0]])
    dm.close()

    path = tmp_path / "no_spectra.zarr"
    _copy_store(base, path)
    root = zarr.open_group(str(path), mode="a")
    del root["spectra"]
    with pytest.raises(ValueError, match="stats/tic requires a spectra/ group"):
        MSIZarrReader(path)

    path2 = tmp_path / "no_ion.zarr"
    _copy_store(base, path2)
    root2 = zarr.open_group(str(path2), mode="a")
    del root2["ion_image"]
    with pytest.raises(
        ValueError, match="stats/mean_spectrum requires an ion_image/ group"
    ):
        MSIZarrReader(path2)


def test_region_mask_output_stats_and_metadata(tmp_path):
    """Masked output recomputes mean_spectrum and updates metadata consistently."""
    source = IonImageDataManagerZarr(
        str(tmp_path / "source.zarr"),
        mode="w",
        width=2,
        height=1,
        temp_dir=str(tmp_path),
        intensity_dtype=np.float32,
        ion_block=16,
    )
    try:
        source.write_flat_batch(
            intensity_flat=np.array([1.0, 2.0, 10.0, 20.0], dtype=np.float32),
            lengths=np.array([2, 2], dtype=np.uint32),
            coordinates=np.array([[1, 1, 1], [2, 1, 1]], dtype=np.uint32),
            mz_axis=np.array([100.0, 200.0]),
        )
        source.finalize()
    finally:
        source.close()

    output: IonImageDataManagerZarr | None = None
    try:
        mask = np.array([[False, True]], dtype=bool)
        output = cast(
            IonImageDataManagerZarr,
            apply_region_mask(source, mask, tmp_path / "masked.zarr"),
        )
        reader = MSIZarrReader(output.filepath)
        try:
            assert reader.root["metadata"].attrs["spectrum_count_num"] == 1
            np.testing.assert_array_equal(
                reader.root["metadata"]["mask"][:], mask
            )
            mean = reader.read_mean_spectrum()
            assert mean is not None
            # 仅保留 (2,1,1) 像素，其谱为 [2, 20]
            np.testing.assert_allclose(mean, [2.0, 20.0])
            reader.validate_full()
        finally:
            reader.close()
    finally:
        if output is not None:
            output.close()
