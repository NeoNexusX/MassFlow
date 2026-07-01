"""Tests for exporting MassFlow ion-image backends to imzML."""

import time
from pathlib import Path

import pytest
import numpy as np

from massflow.data_manager import MSDataManagerImzML, IonImageDataManagerZarr
from massflow.tools.export_imzml import export_imzml
from massflow.tools.logger import get_logger

logger = get_logger("massflow.test_export_imzml")

FILEPATH = Path("data/example.zarr")
ROUND = 5


def _ion_zarr_to_imzml(
    ion_zarr_dm: IonImageDataManagerZarr,
) -> None:
    export_imzml(ion_zarr_dm, "./temp/test_export_imzml.imzML")


def _validation_correct_of_data(
    ion_zarr_dm: IonImageDataManagerZarr,
    imzml_dm: MSDataManagerImzML,
) -> None:
    # Validate that the exported imzML data matches the source Zarr data. include metadata.
    imzml_dm.load_head_data()

    # --- Metadata ----------------------------------------------------------
    zarr_meta = ion_zarr_dm.meta
    imzml_meta = imzml_dm.ms.meta
    assert imzml_meta.continuous is True
    assert imzml_meta.max_count_of_pixels_x == zarr_meta.max_count_of_pixels_x
    assert imzml_meta.max_count_of_pixels_y == zarr_meta.max_count_of_pixels_y

    parser = imzml_dm.parser
    assert ion_zarr_dm.n_pixels == len(parser.coordinates)

    # --- Shared m/z axis ---------------------------------------------------
    zarr_mz = np.asarray(ion_zarr_dm.mz_axis)
    imzml_mz, _ = parser.getspectrum(0)
    imzml_mz = np.asarray(imzml_mz)
    assert ion_zarr_dm.n_ions == zarr_mz.size == imzml_mz.size
    np.testing.assert_allclose(zarr_mz, imzml_mz)

    # --- Coordinate set ----------------------------------------------------
    zarr_coords = np.asarray(ion_zarr_dm.coordinates, dtype=np.int64)
    imzml_coords = np.asarray(parser.coordinates, dtype=np.int64)
    index_by_coord = {tuple(c.tolist()): i for i, c in enumerate(imzml_coords)}
    assert len(index_by_coord) == zarr_coords.shape[0], "imzML coordinates contain duplicates."

    # --- Per-pixel intensities --------------------------------------------
    ion_block = ion_zarr_dm.get_ion_block(0, ion_zarr_dm.n_ions)
    pixel_major = np.asarray(ion_block).T  # (n_pixels, n_ions)
    for pixel_idx, coord in enumerate(zarr_coords):
        imzml_idx = index_by_coord[tuple(coord.tolist())]
        _, intensity = parser.getspectrum(imzml_idx)
        np.testing.assert_allclose(pixel_major[pixel_idx], np.asarray(intensity))


class TestExportImzML:
    """Tests for exporting MassFlow ion-image backends to imzML."""

    @pytest.fixture(scope="class")
    def ion_zarr_dm(self) -> IonImageDataManagerZarr:
        """Fixture that provides an IonImageDataManagerZarr instance for testing."""
        if not FILEPATH.exists():
            pytest.skip(f"Missing test data: {FILEPATH}")

        ion_zarr_dm = IonImageDataManagerZarr(
            filepath=FILEPATH,
            mode="r",
        )
        return ion_zarr_dm

    @pytest.mark.benchmark(timer=time.perf_counter)
    def test_ion_zarr_to_imzml_benchmark(self, benchmark, ion_zarr_dm):
        """Benchmark the export of IonImageDataManagerZarr to imzML."""
        benchmark.pedantic(
            _ion_zarr_to_imzml,
            args=(ion_zarr_dm,),
            rounds=ROUND,
            iterations=1,
            warmup_rounds=0,
        )

    def test_data_validation(self, ion_zarr_dm):
        """Validate that the exported imzML data matches the source Zarr data."""
        imzml_path = Path("./temp/test_export_imzml.imzML")
        if not imzml_path.exists():
            pytest.skip(f"Missing exported imzML file: {imzml_path}")

        imzml_dm = MSDataManagerImzML(filepath=imzml_path)
        _validation_correct_of_data(ion_zarr_dm, imzml_dm)
