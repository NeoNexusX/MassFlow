"""Regression coverage for Binary mask output paths."""

from pathlib import Path

import numpy as np
import pytest

from massflow.data_manager import IonImageDataManagerBinary
from massflow.segmentation.region_mask import apply_region_mask


@pytest.fixture
def binary_source(tmp_path):
    with IonImageDataManagerBinary(
        tmp_path / "source.bin", mode="w", width=2, height=2, temp_dir=tmp_path
    ) as source:
        source.write_flat_batch(
            intensity_flat=np.array([1, 2, 3, 4, 5, 6, 7, 8], dtype=np.float32),
            lengths=np.array([2, 2, 2, 2]),
            coordinates=np.array([[1, 1], [2, 1], [1, 2], [2, 2]]),
            mz_axis=np.array([100.0, 200.0]),
        )
        source.finalize()
        yield source


@pytest.mark.parametrize("path_kind", ["same", "relative", "symlink", "hardlink"])
def test_binary_mask_rejects_source_alias(
    binary_source, tmp_path, monkeypatch, path_kind
):
    source_path = binary_source.filepath
    original = source_path.read_bytes()
    if path_kind == "same":
        output_path = source_path
    elif path_kind == "relative":
        monkeypatch.chdir(tmp_path)
        output_path = Path("source.bin")
    else:
        output_path = tmp_path / "alias.bin"
        if path_kind == "symlink":
            output_path.symlink_to(source_path)
        else:
            output_path.hardlink_to(source_path)

    with pytest.raises(ValueError, match="must not overwrite the source"):
        apply_region_mask(binary_source, np.ones((2, 2), dtype=bool), output_path)

    assert source_path.read_bytes() == original
    np.testing.assert_array_equal(
        binary_source.get_ion_block(0, 2), [[1, 3, 5, 7], [2, 4, 6, 8]]
    )


def test_binary_mask_writes_distinct_output(binary_source, tmp_path):
    original = binary_source.filepath.read_bytes()
    with apply_region_mask(
        binary_source,
        np.array([[True, False], [False, True]]),
        tmp_path / "masked.bin",
    ) as output:
        np.testing.assert_array_equal(output.coordinates, [[1, 1], [2, 2]])
        np.testing.assert_array_equal(output.get_ion_block(0, 2), [[1, 7], [2, 8]])
    assert binary_source.filepath.read_bytes() == original
