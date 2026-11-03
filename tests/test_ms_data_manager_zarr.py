"""Correctness tests for pixel-major MassFlow MSI Zarr stores."""

from __future__ import annotations

from typing import cast

import numpy as np
import pytest
import zarr

from massflow.data_manager import MSDataManagerZarr
from massflow.msi_zarr import (
    MSIZarrReader,
    MSIZarrWriter,
    ZARR_MODE_AUTO,
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
        np.array([100.0, 200.0]),
        np.array([1.0, 10.0], dtype=np.float32),
        np.array([2, 1, 0], dtype=np.uint32),
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


def test_continuous_pixel_major_computes_tic(tmp_path):
    path = tmp_path / "tic_continuous.zarr"
    spectra = [
        (np.array([1.0, 2.0, 3.0], dtype=np.float32), np.array([1, 1, 1], dtype=np.uint32)),
        (np.array([4.0, 5.0, 6.0], dtype=np.float32), np.array([2, 1, 1], dtype=np.uint32)),
    ]
    mz = np.array([100.0, 200.0, 300.0])
    writer = MSIZarrWriter(
        path,
        row_axis="pixel",
        encoding=ZARR_MODE_CONTINUOUS,
        metadata=None,
        width=2,
        height=1,
        compute_tic=True,
        compute_rms=True,
    )
    for intensity, coord in spectra:
        writer.add_spectrum(mz, intensity, coord)
    writer.finalize()
    writer.close()

    expected = np.array([s.sum(dtype=np.float64) for s, _ in spectra], dtype=np.float64)
    expected_rms = np.array(
        [
            np.sqrt(np.mean(np.square(s.astype(np.float64))))
            for s, _ in spectra
        ],
        dtype=np.float64,
    )
    reader = MSIZarrReader(path)
    try:
        tic = reader.read_tic()
        assert tic is not None
        assert tic.dtype == np.float64
        np.testing.assert_allclose(tic, expected, rtol=1e-12)
        tic_array = cast(zarr.Array, reader.root["stats"]["tic"]) # type: ignore
        assert tic_array.shape == (len(spectra),)
        assert tic_array.attrs["source"] == "/spectra"
        assert tic_array.attrs["reduced_axis"] == "mz"
        rms = reader.read_rms()
        assert rms is not None
        assert rms.dtype == np.float64
        np.testing.assert_allclose(rms, expected_rms, rtol=1e-12)
        rms_array = cast(zarr.Array, reader.root["stats"]["rms"]) # type: ignore
        assert rms_array.shape == (len(spectra),)
        assert rms_array.attrs["source"] == "/spectra"
        assert rms_array.attrs["reduced_axis"] == "mz"
        assert rms_array.attrs["denominator"] == "spectrum_points"
        assert reader.n_pixels == len(spectra)
        reader.validate_full()
    finally:
        reader.close()


def test_processed_pixel_major_computes_tic(tmp_path):
    path = tmp_path / "tic_processed.zarr"
    writer = MSIZarrWriter(
        path,
        row_axis="pixel",
        encoding=ZARR_MODE_PROCESSED,
        metadata=None,
        width=3,
        height=1,
        compute_tic=True,
        compute_rms=True,
    )
    writer.write_flat_batch(
        intensity_flat=np.array([1.0, 2.0, 3.0, 4.0], dtype=np.float32),
        lengths=np.array([1, 3], dtype=np.uint32),
        coordinates=np.array([[1, 1, 1], [2, 1, 1]], dtype=np.uint32),
        mz_axis=np.array([100.0, 200.0, 201.0, 300.0]),
    )
    writer.finalize()
    writer.close()

    expected = np.array([1.0, 2.0 + 3.0 + 4.0], dtype=np.float64)
    expected_rms = np.array(
        [1.0, np.sqrt((2.0**2 + 3.0**2 + 4.0**2) / 3.0)],
        dtype=np.float64,
    )
    reader = MSIZarrReader(path)
    try:
        tic = reader.read_tic()
        assert tic is not None
        np.testing.assert_allclose(tic, expected, rtol=1e-12)
        rms = reader.read_rms()
        assert rms is not None
        np.testing.assert_allclose(rms, expected_rms, rtol=1e-12)
        reader.validate_full()
    finally:
        reader.close()


def test_tic_multibatch_matches_row_order(tmp_path):
    path = tmp_path / "tic_multibatch.zarr"
    rng = np.random.default_rng(0)
    n_spectra = 10
    n_mz = 5
    intensities = rng.random((n_spectra, n_mz), dtype=np.float32)
    coords = np.array([[i + 1, 1, 1] for i in range(n_spectra)], dtype=np.uint32)
    mz = np.arange(1, n_mz + 1, dtype=np.float64)
    writer = MSIZarrWriter(
        path,
        row_axis="pixel",
        encoding=ZARR_MODE_CONTINUOUS,
        metadata=None,
        width=n_spectra,
        height=1,
        compute_tic=True,
        compute_rms=True,
        add_buffer_size=3,
    )
    for i in range(n_spectra):
        writer.add_spectrum(mz, intensities[i], coords[i])
    writer.finalize()
    writer.close()

    expected = intensities.sum(axis=1, dtype=np.float64)
    expected_rms = np.sqrt(
        np.mean(np.square(intensities.astype(np.float64)), axis=1)
    )
    reader = MSIZarrReader(path)
    try:
        tic = reader.read_tic()
        assert tic is not None
        assert tic.shape == (n_spectra,)
        np.testing.assert_allclose(tic, expected, rtol=1e-6)
        rms = reader.read_rms()
        assert rms is not None
        assert rms.shape == (n_spectra,)
        np.testing.assert_allclose(rms, expected_rms, rtol=1e-6)
        reader.validate_full()
    finally:
        reader.close()


def test_tic_absent_when_disabled(tmp_path):
    path = tmp_path / "no_tic.zarr"
    writer = MSIZarrWriter(
        path,
        row_axis="pixel",
        encoding=ZARR_MODE_CONTINUOUS,
        metadata=None,
        width=2,
        height=1,
    )
    writer.add_spectrum(
        np.array([100.0, 200.0]),
        np.array([1.0, 2.0], dtype=np.float32),
        np.array([1, 1, 1], dtype=np.uint32),
    )
    writer.finalize()
    writer.close()

    reader = MSIZarrReader(path)
    try:
        assert reader.read_tic() is None
        assert reader.read_rms() is None
        assert "stats" not in reader.root
        reader.validate_full()
    finally:
        reader.close()


def test_compute_tic_rejected_for_ion_major(tmp_path):
    with pytest.raises(ValueError, match="compute_tic is only supported for pixel-major"):
        MSIZarrWriter(
            tmp_path / "ion_tic.zarr",
            row_axis="ion",
            encoding=ZARR_MODE_CONTINUOUS,
            metadata=None,
            width=2,
            height=1,
            compute_tic=True,
        )


def test_compute_rms_rejected_for_ion_major(tmp_path):
    with pytest.raises(ValueError, match="compute_rms is only supported for pixel-major"):
        MSIZarrWriter(
            tmp_path / "ion_rms.zarr",
            row_axis="ion",
            encoding=ZARR_MODE_CONTINUOUS,
            metadata=None,
            width=2,
            height=1,
            compute_rms=True,
        )


# --- 空谱（长度 0）支持：BUG-07 回归 -----------------------------------------


def test_processed_writer_accepts_empty_spectra_roundtrip(tmp_path):
    """首/中/尾空行的 processed 批可写出，读回空数组与原数据一致。"""
    path = tmp_path / "empty_rows.zarr"
    writer = MSIZarrWriter(
        path,
        row_axis="pixel",
        encoding=ZARR_MODE_PROCESSED,
        metadata=None,
        width=4,
        height=1,
    )
    writer.write_flat_batch(
        intensity_flat=np.array([9.0, 8.0, 7.0], dtype=np.float32),
        lengths=np.array([0, 2, 0, 1], dtype=np.uint32),
        coordinates=np.array(
            [[1, 1, 1], [2, 1, 1], [3, 1, 1], [4, 1, 1]], dtype=np.uint32
        ),
        mz_axis=np.array([100.0, 200.0, 300.0]),
    )
    writer.finalize()
    writer.close()

    reader = MSIZarrReader(path)
    try:
        assert reader.encoding == ZARR_MODE_PROCESSED
        offsets = np.asarray(reader.offsets[:], dtype=np.int64)
        np.testing.assert_array_equal(offsets, [0, 0, 2, 2, 3])
        for empty_index in (0, 2):
            mz, intensity = reader.read_spectrum(empty_index)
            assert mz.size == 0 and intensity.size == 0
        mz, intensity = reader.read_spectrum(1)
        np.testing.assert_array_equal(mz, [100.0, 200.0])
        np.testing.assert_array_equal(intensity, [9.0, 8.0])
        mz, intensity = reader.read_spectrum(3)
        np.testing.assert_array_equal(mz, [300.0])
        np.testing.assert_array_equal(intensity, [7.0])
        reader.validate_full()
    finally:
        reader.close()


def test_all_empty_batch_roundtrip_and_stats(tmp_path):
    """全空批：offsets 正确、统计为 0，后续正常批不受影响。"""
    path = tmp_path / "all_empty.zarr"
    writer = MSIZarrWriter(
        path,
        row_axis="pixel",
        encoding=ZARR_MODE_AUTO,
        metadata=None,
        width=3,
        height=1,
        compute_tic=True,
        compute_rms=True,
    )
    writer.write_flat_batch(
        intensity_flat=np.array([], dtype=np.float32),
        lengths=np.array([0, 0], dtype=np.uint32),
        coordinates=np.array([[1, 1, 1], [2, 1, 1]], dtype=np.uint32),
        mz_axis=np.array([], dtype=np.float64),
    )
    writer.write_flat_batch(
        intensity_flat=np.array([1.0, 3.0], dtype=np.float32),
        lengths=np.array([2], dtype=np.uint32),
        coordinates=np.array([[3, 1, 1]], dtype=np.uint32),
        mz_axis=np.array([100.0, 200.0]),
    )
    writer.finalize()
    writer.close()

    reader = MSIZarrReader(path)
    try:
        assert reader.encoding == ZARR_MODE_PROCESSED
        offsets = np.diff(np.asarray(reader.offsets[:], dtype=np.int64))
        np.testing.assert_array_equal(offsets, [0, 0, 2])
        np.testing.assert_array_equal(reader.read_tic(), [0.0, 0.0, 4.0])
        np.testing.assert_allclose(
            reader.read_rms(), [0.0, 0.0, np.sqrt(5.0)], rtol=1e-12
        )
        reader.validate_full()
    finally:
        reader.close()


def test_tic_rms_zero_for_empty_rows_in_any_position(tmp_path):
    """空行在首/中/尾时 TIC=RMS=0，非空行与手算一致（尾空行曾致 reduceat 越界）。"""
    path = tmp_path / "empty_stats.zarr"
    writer = MSIZarrWriter(
        path,
        row_axis="pixel",
        encoding=ZARR_MODE_PROCESSED,
        metadata=None,
        width=5,
        height=1,
        compute_tic=True,
        compute_rms=True,
    )
    values = [1.0, 2.0, 3.0, 4.0]
    writer.write_flat_batch(
        intensity_flat=np.array(values, dtype=np.float32),
        lengths=np.array([0, 2, 0, 2, 0], dtype=np.uint32),
        coordinates=np.array(
            [[i + 1, 1, 1] for i in range(5)], dtype=np.uint32
        ),
        mz_axis=np.array([100.0, 200.0, 300.0, 400.0]),
    )
    writer.finalize()
    writer.close()

    reader = MSIZarrReader(path)
    try:
        tic = reader.read_tic()
        rms = reader.read_rms()
        np.testing.assert_allclose(
            tic, [0.0, 3.0, 0.0, 7.0, 0.0], rtol=1e-12
        )
        np.testing.assert_allclose(
            rms,
            [0.0, np.sqrt(2.5), 0.0, np.sqrt(12.5), 0.0],
            rtol=1e-12,
        )
        reader.validate_full()
    finally:
        reader.close()


@pytest.mark.parametrize("leading_intensity", [1e16, np.nan, np.inf])
@pytest.mark.parametrize(
    ("compute_tic", "compute_rms"), [(True, False), (False, True), (True, True)]
)
def test_spectrum_stats_do_not_leak_between_rows(
    tmp_path, leading_intensity, compute_tic, compute_rms
):
    """Large or nonfinite values must not corrupt later spectra or empty rows."""
    path = tmp_path / "independent_stats.zarr"
    writer = MSIZarrWriter(
        path,
        row_axis="pixel",
        encoding=ZARR_MODE_PROCESSED,
        metadata=None,
        width=5,
        height=1,
        intensity_dtype=np.float64,
        compute_tic=compute_tic,
        compute_rms=compute_rms,
    )
    try:
        writer.write_flat_batch(
            intensity_flat=np.array([leading_intensity, 1.0], dtype=np.float64),
            lengths=np.array([0, 1, 0, 1, 0], dtype=np.int64),
            coordinates=np.array([[i + 1, 1, 1] for i in range(5)]),
            mz_axis=np.array([100.0, 200.0]),
        )
        writer.finalize()
    finally:
        writer.close()

    expected = np.array([0.0, leading_intensity, 0.0, 1.0, 0.0])
    reader = MSIZarrReader(path)
    try:
        if compute_tic:
            np.testing.assert_allclose(reader.read_tic(), expected, equal_nan=True)
        if compute_rms:
            np.testing.assert_allclose(reader.read_rms(), expected, equal_nan=True)
    finally:
        reader.close()


def test_auto_inference_prefers_processed_with_empty_rows(tmp_path):
    """lengths=[4,0] 且 mz 恰为 4 点：歧义数据按唯一合法解释定为 processed。"""
    path = tmp_path / "ambiguous.zarr"
    writer = MSIZarrWriter(
        path,
        row_axis="pixel",
        encoding=ZARR_MODE_AUTO,
        metadata=None,
        width=2,
        height=1,
    )
    mz = np.array([100.0, 150.0, 200.0, 250.0])
    writer.write_flat_batch(
        intensity_flat=np.array([1.0, 2.0, 3.0, 4.0], dtype=np.float32),
        lengths=np.array([4, 0], dtype=np.uint32),
        coordinates=np.array([[1, 1, 1], [2, 1, 1]], dtype=np.uint32),
        mz_axis=mz,
    )
    writer.finalize()
    writer.close()

    reader = MSIZarrReader(path)
    try:
        assert reader.encoding == ZARR_MODE_PROCESSED
        mz0, intensity0 = reader.read_spectrum(0)
        np.testing.assert_array_equal(mz0, mz)
        np.testing.assert_array_equal(intensity0, [1.0, 2.0, 3.0, 4.0])
        mz1, intensity1 = reader.read_spectrum(1)
        assert mz1.size == 0 and intensity1.size == 0
        reader.validate_full()
    finally:
        reader.close()


def test_auto_inference_shared_axis_with_empty_reports_contract_error(tmp_path):
    """共享轴批混入空行：推断为 continuous 后被等长契约明确拒绝。"""
    writer = MSIZarrWriter(
        tmp_path / "mixed_shared.zarr",
        row_axis="pixel",
        encoding=ZARR_MODE_AUTO,
        metadata=None,
        width=3,
        height=1,
    )
    with pytest.raises(ValueError, match="must match axes/mz length"):
        writer.write_flat_batch(
            intensity_flat=np.array([1.0, 2.0, 3.0, 4.0], dtype=np.float32),
            lengths=np.array([2, 0, 2], dtype=np.uint32),
            coordinates=np.array(
                [[1, 1, 1], [2, 1, 1], [3, 1, 1]], dtype=np.uint32
            ),
            mz_axis=np.array([100.0, 200.0]),
        )


def test_negative_lengths_rejected(tmp_path):
    writer = MSIZarrWriter(
        tmp_path / "negative.zarr",
        row_axis="pixel",
        encoding=ZARR_MODE_PROCESSED,
        metadata=None,
        width=2,
        height=1,
    )
    with pytest.raises(ValueError, match="non-negative"):
        writer.write_flat_batch(
            intensity_flat=np.array([1.0, 2.0], dtype=np.float32),
            lengths=np.array([-1, 3], dtype=np.int64),
            coordinates=np.array([[1, 1, 1], [2, 1, 1]], dtype=np.uint32),
            mz_axis=np.array([100.0, 200.0]),
        )


def test_add_spectrum_empty_roundtrip(tmp_path):
    """add_spectrum 传入空谱：flush 后按 processed 写出并可读回。"""
    path = tmp_path / "add_empty.zarr"
    writer = MSIZarrWriter(
        path,
        row_axis="pixel",
        encoding=ZARR_MODE_AUTO,
        metadata=None,
        width=2,
        height=1,
    )
    writer.add_spectrum(
        np.array([], dtype=np.float64),
        np.array([], dtype=np.float32),
        np.array([1, 1, 1], dtype=np.uint32),
    )
    writer.add_spectrum(
        np.array([100.0, 200.0]),
        np.array([5.0, 6.0], dtype=np.float32),
        np.array([2, 1, 1], dtype=np.uint32),
    )
    writer.finalize()
    writer.close()

    reader = MSIZarrReader(path)
    try:
        assert reader.encoding == ZARR_MODE_PROCESSED
        mz, intensity = reader.read_spectrum(0)
        assert mz.size == 0 and intensity.size == 0
        mz, intensity = reader.read_spectrum(1)
        np.testing.assert_array_equal(mz, [100.0, 200.0])
        np.testing.assert_array_equal(intensity, [5.0, 6.0])
        reader.validate_full()
    finally:
        reader.close()


def test_pipeline_spectra_zarr_with_background_pixels(tmp_path):
    """BUG-07 端到端：峰提取后含空谱，spectra_zarr 输出不再崩溃。"""
    from massflow.data_manager import MSDataManagerImzML
    from massflow.preprocess import Preprocessor
    from massflow.tools.stream_imzml_writer import ImzMLWriter

    rng = np.random.default_rng(0)
    mz_axis = np.linspace(100.0, 300.0, 400)
    imzml_path = tmp_path / "background.imzML"
    writer = ImzMLWriter(
        str(imzml_path), mode="continuous", spec_type="profile"
    )
    for i in range(32):
        background = i % 4 == 0
        intensity = (
            np.zeros(400, np.float32)
            if background
            else (
                rng.random(400)
                + 50 * np.exp(-((mz_axis - 200.0) / 2.0) ** 2)
            ).astype(np.float32)
        )
        writer.add_spectrum(
            mz_axis, intensity, (i % 8 + 1, i // 8 + 1, 1)
        )
    writer.close()

    source = MSDataManagerImzML(filepath=str(imzml_path))
    source.load_head_data()
    output = Preprocessor(
        source,
        temp_dir=str(tmp_path / "temp"),
        batch_size=8,
        output_format="spectra_zarr",
        output_path=str(tmp_path / "out.zarr"),
    ).peak_pick().start()

    reader = MSIZarrReader(tmp_path / "out.zarr")
    try:
        assert reader.encoding == ZARR_MODE_PROCESSED
        offsets = np.diff(np.asarray(reader.offsets[:], dtype=np.int64))
        assert int((offsets == 0).sum()) == 8  # 每 4 个像素一个背景
        assert offsets.size == 32
        mz, intensity = reader.read_spectrum(0)
        assert mz.size == 0 and intensity.size == 0
        mz, intensity = reader.read_spectrum(1)
        assert mz.size > 0 and intensity.size == mz.size
        reader.validate_full()
    finally:
        reader.close()
        output.close()
        source.close()
