from pathlib import Path

import numpy as np
import pytest

import matplotlib

matplotlib.use("Agg")

from massflow.data_manager import (
    IonImageDataManagerBinary,
    IonImageDataManagerH5,
    IonImageDataManagerZarr,
    MSDataManagerImzML,
)
from massflow.tools.plot_ion_image import (
    plot_ion_image,
    plot_tic_image,
    read_ion_image,
    read_tic_image,
)

CONTINUOUS_IMZML = Path("./data/example_continuous.imzML")
PROCESSED_IMZML = Path("./data/example_processed.imzML")
CONTINUOUS_TARGET_LOCS = ([1, 39], [2, 40])
PROCESSED_TARGET_LOCS = ([1, 1], [2, 1])


def _require_imzml_pair(path: Path) -> None:
    if not path.exists() or not path.with_suffix(".ibd").exists():
        pytest.skip(f"Missing real imzML test data: {path} and matching .ibd")


def _image_shape(dm: MSDataManagerImzML) -> tuple[int, int]:
    return int(dm.ms.meta.max_count_of_pixels_y), int(dm.ms.meta.max_count_of_pixels_x) # type: ignore


def _scatter_values(coordinates: np.ndarray, values: np.ndarray, shape: tuple[int, int]) -> np.ndarray:
    image = np.zeros(shape, dtype=values.dtype)
    coords = np.asarray(coordinates)
    xs = coords[:, 0].astype(np.int64, copy=False)
    ys = coords[:, 1].astype(np.int64, copy=False)

    zero_based = bool((xs.size > 0 and int(xs.min()) == 0) or (ys.size > 0 and int(ys.min()) == 0))
    if not zero_based:
        xs = xs - 1
        ys = ys - 1

    valid = (xs >= 0) & (xs < shape[1]) & (ys >= 0) & (ys < shape[0])
    image[ys[valid], xs[valid]] = values[valid]
    return image


def _sum_by_lengths(intensity_flat: np.ndarray, lengths: np.ndarray) -> np.ndarray:
    values = np.zeros(lengths.size, dtype=np.result_type(intensity_flat.dtype, np.float64))
    offset = 0
    for index, length_value in enumerate(lengths):
        length = int(length_value)
        end = offset + length
        values[index] = np.sum(intensity_flat[offset:end])
        offset = end
    return values


def _expected_tic_image(dm: MSDataManagerImzML) -> np.ndarray:
    coord_parts = []
    value_parts = []
    for _mz_axis, intensity_flat, lengths, coordinates in dm.flat_generator(
        batch_size=4,
        include_mz=False,
        max_threads=2,
    ):
        coord_parts.append(np.asarray(coordinates))
        value_parts.append(_sum_by_lengths(np.asarray(intensity_flat), np.asarray(lengths, dtype=np.int64)))

    return _scatter_values(np.vstack(coord_parts), np.concatenate(value_parts), _image_shape(dm))


def _expected_shared_ion_image(dm: MSDataManagerImzML, ion_idx: int) -> np.ndarray:
    coord_parts = []
    value_parts = []
    for mz_axis, intensity_flat, lengths, coordinates in dm.flat_generator(
        batch_size=4,
        include_mz=True,
        max_threads=2,
    ):
        mz_axis = np.asarray(mz_axis)
        lengths = np.asarray(lengths, dtype=np.int64)
        intensity_matrix = np.asarray(intensity_flat).reshape(lengths.size, mz_axis.size)
        coord_parts.append(np.asarray(coordinates))
        value_parts.append(intensity_matrix[:, ion_idx])

    return _scatter_values(np.vstack(coord_parts), np.concatenate(value_parts), _image_shape(dm))


def _expected_processed_ion_image(
    dm: MSDataManagerImzML,
    mz: float,
    mz_tolerance: float | None = None,
) -> np.ndarray:
    coord_parts = []
    value_parts = []
    if mz_tolerance is not None:
        low = mz - mz_tolerance
        high = mz + mz_tolerance
    else:
        low = high = None

    for mz_flat, intensity_flat, lengths, coordinates in dm.flat_generator(
        batch_size=2,
        include_mz=True,
        max_threads=2,
    ):
        mz_flat = np.asarray(mz_flat)
        intensity_flat = np.asarray(intensity_flat)
        lengths = np.asarray(lengths, dtype=np.int64)
        values = np.zeros(lengths.size, dtype=np.result_type(intensity_flat.dtype, np.float64))

        offset = 0
        for index, length_value in enumerate(lengths):
            length = int(length_value)
            end = offset + length
            mz_slice = mz_flat[offset:end]
            intensity_slice = intensity_flat[offset:end]
            if length > 0:
                if mz_tolerance is None:
                    nearest = int(np.nanargmin(np.abs(mz_slice.astype(float, copy=False) - float(mz))))
                    values[index] = intensity_slice[nearest]
                else:
                    mask = (mz_slice >= low) & (mz_slice <= high)
                    values[index] = np.sum(intensity_slice[mask])
            offset = end

        coord_parts.append(np.asarray(coordinates))
        value_parts.append(values)

    return _scatter_values(np.vstack(coord_parts), np.concatenate(value_parts), _image_shape(dm))


def _write_backend_from_real_subset(
    dm: MSDataManagerImzML,
    manager_cls,
    output_path: Path,
    tmp_path: Path,
    *,
    ion_count: int = 8,
) -> np.ndarray:
    shape = _image_shape(dm)
    expected_stack = np.zeros((ion_count, shape[0], shape[1]), dtype=np.float64)

    with manager_cls(
        filepath=output_path,
        mode="w",
        width=shape[1],
        height=shape[0],
        temp_dir=tmp_path,
        mz_dtype=dm.mz_dtype,
        intensity_dtype=dm.intensity_dtype,
    ) as backend_dm:
        for mz_axis, intensity_flat, lengths, coordinates in dm.flat_generator(
            batch_size=4,
            include_mz=True,
            max_threads=2,
        ):
            mz_axis = np.asarray(mz_axis)
            lengths = np.asarray(lengths, dtype=np.int64)
            coordinates = np.asarray(coordinates)
            n_pixels = lengths.size
            intensity_matrix = np.asarray(intensity_flat).reshape(n_pixels, mz_axis.size)
            trimmed_matrix = intensity_matrix[:, :ion_count]
            trimmed_lengths = np.full(n_pixels, ion_count, dtype=np.uint32)

            backend_dm.write_flat_batch(
                intensity_flat=trimmed_matrix.reshape(-1),
                lengths=trimmed_lengths,
                coordinates=coordinates,
                mz_axis=mz_axis[:ion_count],
            )

            for pixel_idx, coord in enumerate(coordinates):
                x = int(coord[0]) - 1
                y = int(coord[1]) - 1
                expected_stack[:, y, x] = trimmed_matrix[pixel_idx]

        backend_dm.finalize()

    return expected_stack


class TestPlotIonImage:
    """
    Test suite for ion image plotting functionality.
    """
    @pytest.fixture(scope="class")
    def continuous_dm(self):
        _require_imzml_pair(CONTINUOUS_IMZML)
        dm = MSDataManagerImzML(filepath=str(CONTINUOUS_IMZML), target_locs=CONTINUOUS_TARGET_LOCS)
        dm.load_head_data()
        if len(dm.ms) == 0:
            dm.close()
            pytest.skip(f"No spectra found in target_locs={CONTINUOUS_TARGET_LOCS}")
        yield dm
        dm.close()

    @pytest.fixture(scope="class")
    def processed_dm(self):
        _require_imzml_pair(PROCESSED_IMZML)
        dm = MSDataManagerImzML(filepath=str(PROCESSED_IMZML), target_locs=PROCESSED_TARGET_LOCS)
        dm.load_head_data()
        if len(dm.ms) == 0:
            dm.close()
            pytest.skip(f"No spectra found in target_locs={PROCESSED_TARGET_LOCS}")
        yield dm
        dm.close()

    def test_read_tic_image_from_real_imzml_matches_flat_generator(self, continuous_dm):
        image = read_tic_image(continuous_dm, batch_size=4, max_threads=2)
        expected = _expected_tic_image(continuous_dm)

        np.testing.assert_allclose(image, expected, rtol=1e-6, atol=1e-6)
        assert image.shape == _image_shape(continuous_dm)
        assert np.nanmax(image) > 0

    def test_read_ion_image_from_real_imzml_maps_mz_to_nearest_shared_axis(self, continuous_dm):
        ion_idx = 5
        mz_axis = np.asarray(continuous_dm.ms.shared_mz_list)
        requested_mz = float(mz_axis[ion_idx] + (mz_axis[ion_idx + 1] - mz_axis[ion_idx]) * 0.1)

        image = read_ion_image(continuous_dm, mz=requested_mz, batch_size=4, max_threads=2)
        expected = _expected_shared_ion_image(continuous_dm, ion_idx)

        np.testing.assert_allclose(image, expected, rtol=1e-6, atol=1e-6)

    def test_read_ion_image_requires_mz(self, continuous_dm):
        with pytest.raises(ValueError, match="'mz' is required"):
            read_ion_image(continuous_dm)

    def test_processed_real_imzml_nearest_without_mz_tolerance_matches_flat_generator(self, processed_dm):
        mz_flat, _intensity_flat, _lengths, _coordinates = next(
            processed_dm.flat_generator(batch_size=2, include_mz=True, max_threads=2)
        )
        requested_mz = float(np.asarray(mz_flat)[0])

        image = read_ion_image(
            processed_dm,
            mz=requested_mz,
            batch_size=2,
            max_threads=2,
        )
        expected = _expected_processed_ion_image(processed_dm, requested_mz)

        np.testing.assert_allclose(image, expected, rtol=1e-6, atol=1e-6)

    def test_processed_real_imzml_mz_window_matches_flat_generator(self, processed_dm):
        mz_flat, _intensity_flat, _lengths, _coordinates = next(
            processed_dm.flat_generator(batch_size=2, include_mz=True, max_threads=2)
        )
        requested_mz = float(np.asarray(mz_flat)[0])
        mz_tolerance = 0.01

        image = read_ion_image(
            processed_dm,
            mz=requested_mz,
            mz_tolerance=mz_tolerance,
            batch_size=2,
            max_threads=2,
        )
        expected = _expected_processed_ion_image(processed_dm, requested_mz, mz_tolerance)

        np.testing.assert_allclose(image, expected, rtol=1e-6, atol=1e-6)

    def test_structure_only_tic_plot_from_real_imzml_only_changes_plot_style(self, continuous_dm, tmp_path):
        plot_path = tmp_path / "tic_image.png"

        plot_tic_image(
            continuous_dm,
            structure_only=True,
            output_dir=tmp_path,
            show=False,
            dpi=123,
        )

        assert plot_path.exists()

    @pytest.mark.parametrize(
        ("manager_cls", "suffix"),
        [
            (IonImageDataManagerBinary, ".iib"),
            (IonImageDataManagerH5, ".h5"),
            (IonImageDataManagerZarr, ".zarr"),
        ],
    )
    def test_real_derived_ion_image_backend_managers_read_and_plot(
        self,
        continuous_dm,
        tmp_path,
        manager_cls,
        suffix,
    ):
        output_path = tmp_path / f"real_subset{suffix}"
        expected_stack = _write_backend_from_real_subset(
            continuous_dm,
            manager_cls,
            output_path,
            tmp_path,
        )

        mz_axis = np.asarray(continuous_dm.ms.shared_mz_list)
        with manager_cls(filepath=output_path, mode="r") as backend_dm:
            image = read_ion_image(backend_dm, mz=float(mz_axis[1]))
            tic_image = read_tic_image(backend_dm)
            plot_path = tmp_path / f"plot{suffix}.png"
            plot_ion_image(backend_dm, mz=float(mz_axis[1]), save_path=plot_path, show=False, dpi=123)
            tic_cover_dir = tmp_path / f"tic_cover{suffix}"
            plot_tic_image(backend_dm, structure_only=True, output_dir=tic_cover_dir, show=False, dpi=123)
            tic_cover_path = tic_cover_dir / "tic_image.png"

        np.testing.assert_allclose(image, expected_stack[1], rtol=1e-5, atol=1e-5)
        np.testing.assert_allclose(tic_image, np.sum(expected_stack, axis=0), rtol=1e-5, atol=1e-5)
        assert plot_path.exists()
        assert tic_cover_path.exists()
