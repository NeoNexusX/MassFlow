"""Regression coverage for H5 coordinate bases and bounds."""

import numpy as np
import pytest

from massflow.data_manager import IonImageDataManagerH5


@pytest.mark.parametrize("infer_dimensions", [False, True])
@pytest.mark.parametrize(
    "coordinates, inferred_shape, expected_image",
    [
        (
            [[1, 1], [2, 1], [1, 2], [2, 2]],
            (2, 2),
            [[1, 3, 0], [5, 7, 0], [0, 0, 0]],
        ),
        (
            [[1, 0], [2, 0], [1, 1], [2, 1]],
            (2, 3),
            [[0, 1, 3], [0, 5, 7], [0, 0, 0]],
        ),
        (
            [[0, 1], [1, 1], [0, 2], [1, 2]],
            (3, 2),
            [[0, 0, 0], [1, 3, 0], [5, 7, 0]],
        ),
        (
            [[0, 0], [1, 0], [0, 1], [1, 1]],
            (2, 2),
            [[1, 3, 0], [5, 7, 0], [0, 0, 0]],
        ),
    ],
)
def test_h5_coordinate_roundtrip(
    tmp_path, infer_dimensions, coordinates, inferred_shape, expected_image
):
    path = tmp_path / "coordinates.h5"
    height, width = inferred_shape if infer_dimensions else (3, 3)
    expected_image = np.array(expected_image, dtype=np.float32)[:height, :width]
    expected_block = np.array([[1, 3, 5, 7], [2, 4, 6, 8]], dtype=np.float32)
    with IonImageDataManagerH5(
        path,
        mode="w",
        width=None if infer_dimensions else width,
        height=None if infer_dimensions else height,
        temp_dir=tmp_path,
    ) as writer:
        writer.write_flat_batch(
            intensity_flat=np.array([1, 2, 3, 4, 5, 6, 7, 8], dtype=np.float32),
            lengths=np.array([2, 2, 2, 2]),
            coordinates=np.array(coordinates),
            mz_axis=np.array([100.0, 200.0]),
        )
        writer.finalize()
        np.testing.assert_array_equal(writer.get_ion_image(0), expected_image)
        np.testing.assert_array_equal(writer.get_ion_block(0, 2), expected_block)

    with IonImageDataManagerH5(path, mode="r") as reader:
        assert (reader.height, reader.width) == (height, width)
        assert reader.n_pixels == len(coordinates)
        np.testing.assert_array_equal(reader.coordinates, coordinates)
        np.testing.assert_array_equal(reader.get_ion_image(0), expected_image)
        np.testing.assert_array_equal(reader.get_ion_vector(1), expected_block[1])
        np.testing.assert_array_equal(reader.get_ion_block(0, 2), expected_block)


@pytest.mark.parametrize("coordinates", [[[1, 0], [3, 1]], [[1, 1], [2, 3]]])
def test_h5_out_of_bounds_preserves_existing_output(tmp_path, coordinates):
    path = tmp_path / "existing.h5"
    original = b"existing output must not be truncated"
    path.write_bytes(original)
    with IonImageDataManagerH5(
        path, mode="w", width=3, height=2, temp_dir=tmp_path
    ) as writer:
        writer.write_flat_batch(
            intensity_flat=np.array([1, 2], dtype=np.float32),
            lengths=np.array([1, 1]),
            coordinates=np.array(coordinates),
            mz_axis=np.array([100.0]),
        )
        with pytest.raises(ValueError, match="outside the H5 image bounds"):
            writer.finalize()
    assert path.read_bytes() == original


def test_h5_rejects_negative_coordinates_before_unsigned_conversion(tmp_path):
    path = tmp_path / "negative.h5"
    with IonImageDataManagerH5(path, mode="w", temp_dir=tmp_path) as writer:
        with pytest.raises(ValueError, match="coordinates must be non-negative"):
            writer.write_flat_batch(
                intensity_flat=np.array([1], dtype=np.float32),
                lengths=np.array([1]),
                coordinates=np.array([[1, -1]]),
                mz_axis=np.array([100.0]),
            )
    assert not path.exists()
