"""Tests for exporting MassFlow ion-image backends to imzML."""

import hashlib
from pathlib import Path

import numpy as np
import pytest
from pyimzml.ImzMLParser import ImzMLParser

from massflow.data_manager import (
    IonImageDataManagerBinary,
    IonImageDataManagerH5,
    IonImageDataManagerZarr,
)
from massflow.tools.export_imzml import export_imzml


MZ_AXIS = np.array([100.0, 200.0, 300.0], dtype=np.float64)
COORDINATES = np.array([[1, 1], [2, 1], [1, 2]], dtype=np.uint32)
PIXEL_SPECTRA = np.array(
    [
        [1.0, 2.0, 3.0],
        [10.0, 20.0, 30.0],
        [100.0, 200.0, 300.0],
    ],
    dtype=np.float32,
)


def _write_source(backend_class, path: Path, temp_dir: Path) -> None:
    with backend_class(
        path,
        mode="w",
        width=2,
        height=2,
        temp_dir=temp_dir,
        ion_block=16,
    ) as data_manager:
        data_manager.write_flat_batch(
            intensity_flat=PIXEL_SPECTRA.ravel(),
            lengths=np.full(len(COORDINATES), len(MZ_AXIS), dtype=np.uint32),
            coordinates=COORDINATES,
            mz_axis=MZ_AXIS,
        )
        data_manager.finalize()


@pytest.mark.parametrize(
    ("backend_class", "suffix"),
    [
        (IonImageDataManagerBinary, ".iib"),
        (IonImageDataManagerH5, ".h5"),
        (IonImageDataManagerZarr, ".zarr"),
    ],
)
def test_export_imzml_roundtrips_all_backends(tmp_path, backend_class, suffix):
    source_path = tmp_path / f"source{suffix}"
    output_path = tmp_path / "exported.imzML"
    _write_source(backend_class, source_path, tmp_path)

    imzml_path, ibd_path = export_imzml(
        source_path,
        output_path,
        temp_dir=tmp_path,
        ion_block_size=2,
        spectrum_batch_size=2,
    )

    assert imzml_path == output_path
    assert ibd_path == output_path.with_suffix(".ibd")
    assert (
        hashlib.sha1(ibd_path.read_bytes()).hexdigest().upper()
        in imzml_path.read_text() # pylint: disable=W1514
    )
    parser = ImzMLParser(str(imzml_path))
    assert parser.coordinates == [tuple(coord) + (1,) for coord in COORDINATES.tolist()]
    for index in range(len(COORDINATES)):
        mz_axis, intensity = parser.getspectrum(index)
        np.testing.assert_array_equal(mz_axis, MZ_AXIS)
        np.testing.assert_array_equal(intensity, PIXEL_SPECTRA[index])
