"""Tests for exporting MassFlow ion-image backends to imzML."""

# pylint: disable=too-many-arguments,too-many-positional-arguments,too-many-locals

import time
import uuid
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np
import pytest
from pyimzml.ImzMLParser import ImzMLParser
from pyimzml.ImzMLWriter import ImzMLWriter as OfficialImzMLWriter

from massflow.data_manager import (
    IonImageDataManagerBinary,
    IonImageDataManagerH5,
    IonImageDataManagerZarr,
    MSDataManagerImzML,
    MSDataManagerZarr,
)
from massflow.tools.export_imzml import export_imzml
from massflow.tools.logger import get_logger
from massflow.tools.stream_imzml_writer import PLACEHOLDER_SHA1, ImzMLWriter

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


# ---------------------------------------------------------------------------
# Self-contained export fixtures and binary/XML assertions
# ---------------------------------------------------------------------------

_NS = "{http://psi.hupo.org/ms/mzml}"
_COORDINATES = np.array([[2, 2, 1], [3, 2, 1], [4, 2, 1], [2, 3, 1], [3, 3, 1]])
_MZ = np.array([100.0, 200.0, 300.0], dtype=np.float64)
_VALUES = np.arange(1, 16, dtype=np.float32).reshape(5, 3)


def _write_ion_source(tmp_path, backend=IonImageDataManagerZarr):
    suffix = {
        IonImageDataManagerZarr: ".zarr",
        IonImageDataManagerBinary: ".bin",
        IonImageDataManagerH5: ".h5",
    }[backend]
    source = backend(
        tmp_path / ("ion_source" + suffix),
        mode="w",
        width=6,
        height=5,
        temp_dir=tmp_path,
        ion_block=2,
    )
    source.write_flat_batch(
        intensity_flat=_VALUES.reshape(-1),
        lengths=np.full(5, 3, dtype=np.int64),
        coordinates=_COORDINATES,
        mz_axis=_MZ,
    )
    source.finalize()
    return source


def _write_spectra_source(tmp_path, mode):
    source = MSDataManagerZarr(
        filepath=tmp_path / "spectra.zarr",
        mode="w",
        storage_mode=mode,
        temp_dir=tmp_path,
    )
    source.ms.meta.max_count_of_pixels_x = 6
    source.ms.meta.max_count_of_pixels_y = 5
    lengths = (
        np.full(5, 3, dtype=np.int64)
        if mode == "continuous"
        else np.array([2, 4, 1, 3, 2])
    )
    mz = (
        _MZ
        if mode == "continuous"
        else np.arange(lengths.sum(), dtype=np.float64) + 100
    )
    intensity = (
        _VALUES.reshape(-1)
        if mode == "continuous"
        else np.arange(lengths.sum(), dtype=np.float32) + 10
    )
    source.write_flat_batch(
        intensity_flat=intensity,
        lengths=lengths,
        coordinates=_COORDINATES,
        mz_axis=mz,
    )
    source.finalize()
    source.close()
    return MSDataManagerZarr(filepath=tmp_path / "spectra.zarr"), mz, intensity, lengths


def _assert_spectra(path, coordinates, mz, intensity, lengths, indices, mode):
    offsets = np.concatenate(([0], np.cumsum(lengths)))
    with ImzMLParser(str(path)) as parser:
        np.testing.assert_array_equal(parser.coordinates, coordinates[indices])
        for output_index, source_index in enumerate(indices):
            start, end = offsets[source_index : source_index + 2]
            actual_mz, actual_intensity = parser.getspectrum(output_index)
            np.testing.assert_array_equal(
                actual_mz, mz if mode == "continuous" else mz[start:end]
            )
            np.testing.assert_array_equal(actual_intensity, intensity[start:end])
        if mode == "continuous":
            assert len(set(parser.mzOffsets)) == 1
        else:
            for previous, current, length in zip(
                parser.mzOffsets, parser.mzOffsets[1:], parser.mzLengths
            ):
                assert current >= previous + length * np.dtype(np.float64).itemsize
        root = parser.root
        assert root is not None
        checksum = root.find(".//" + _NS + "cvParam[@accession='IMS:1000091']")
        assert checksum is not None
        assert checksum.attrib["value"] == PLACEHOLDER_SHA1
        uuid_param = root.find(".//" + _NS + "cvParam[@accession='IMS:1000080']")
        assert uuid_param is not None
        ibd_uuid = uuid.UUID(uuid_param.attrib["value"])
        with path.with_suffix(".ibd").open("rb") as binary:
            assert binary.read(16) == ibd_uuid.bytes
    # A new tree is needed because the parser prunes spectrum nodes.
    root = ET.parse(path).getroot()
    spectra = root.findall(".//" + _NS + "spectrum")
    assert len(spectra) == len(indices)
    for spectrum in spectra:
        for accession in (
            "MS:1000528",
            "MS:1000527",
            "MS:1000504",
            "MS:1000505",
            "MS:1000285",
        ):
            statistic = spectrum.find(_NS + "cvParam[@accession='" + accession + "']")
            assert statistic is not None
            assert float(statistic.attrib["value"]) == 0


# ---------------------------------------------------------------------------
# Ion-image transpose and spectra stream regressions
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "backend",
    [IonImageDataManagerZarr, IonImageDataManagerBinary, IonImageDataManagerH5],
)
@pytest.mark.parametrize("batch_size", [1, 2])
@pytest.mark.parametrize("masked", [False, True])
def test_ion_export_batches(tmp_path, backend, batch_size, masked, monkeypatch):
    with _write_ion_source(tmp_path, backend) as source:
        mask = None
        indices = np.arange(5)
        if masked:
            mask = np.zeros((5, 6), dtype=bool)
            mask[1, 1] = mask[1, 3] = mask[2, 2] = True
            indices = np.array([0, 2, 4])
        reads = []
        get_block = source.get_ion_block

        def record_read(start, end):
            reads.append((start, end))
            return get_block(start, end)

        monkeypatch.setattr(source, "get_ion_block", record_read)
        temp_dir = tmp_path / "new" / "transpose"
        path = export_imzml(
            source,
            tmp_path / "out.imzML",
            mask=mask,
            batch_size=batch_size,
            temp_dir=temp_dir,
        )
        assert temp_dir.is_dir()
        assert not list(temp_dir.iterdir())
        assert [ion for start, end in reads for ion in range(start, end)] == [0, 1, 2]
        _assert_spectra(
            path,
            _COORDINATES,
            _MZ,
            _VALUES.reshape(-1),
            np.full(5, 3),
            indices,
            "continuous",
        )


@pytest.mark.parametrize("mode", ["continuous", "processed"])
@pytest.mark.parametrize("selection", ["full", "all", "partial"])
def test_spectra_zarr_export(tmp_path, mode, selection):
    source, mz, intensity, lengths = _write_spectra_source(tmp_path, mode)
    with source:
        mask = None if selection == "full" else np.ones((5, 6), dtype=bool)
        indices = np.arange(5)
        if selection == "partial":
            mask = np.zeros((5, 6), dtype=bool)
            mask[1, 1] = mask[1, 3] = True
            indices = np.array([0, 2])  # partial, partial, and empty batches
        path = export_imzml(source, tmp_path / "out.imzML", mask=mask, batch_size=2)
        _assert_spectra(path, _COORDINATES, mz, intensity, lengths, indices, mode)


def test_zarr_to_imzml_copy_meta_cross_type(tmp_path):
    """copy_meta must copy shared fields from a MetaDataZarr source (LOGIC-01)."""
    source, _, _, _ = _write_spectra_source(tmp_path, "processed")
    source.extract_metadata()  # load on-disk metadata, as the pipeline read path does
    source.ms.meta.pixel_size_x = 25.0
    source.ms.meta.pixel_size_y = 25.0
    source.ms.meta.instrument_model = "test instrument"

    with MSDataManagerImzML(filepath=str(tmp_path / "out.imzML")) as dm:
        dm.copy_meta(source)
        assert dm.ms.meta.pixel_size_x == 25.0
        assert dm.ms.meta.pixel_size_y == 25.0
        assert dm.ms.meta.instrument_model == "test instrument"
        assert dm.ms.meta.max_count_of_pixels_x == 6
        # Fields excluded by copy must not leak from the source.
        assert dm.ms.meta.processed is None


def _write_imzml_source(tmp_path, mode, coordinates=_COORDINATES):
    path = tmp_path / "source.imzML"
    lengths = np.full(5, 3) if mode == "continuous" else np.array([2, 4, 1, 3, 2])
    mz = (
        _MZ
        if mode == "continuous"
        else np.arange(lengths.sum(), dtype=np.float64) + 100
    )
    intensity = np.arange(lengths.sum(), dtype=np.float32) + 10
    offsets = np.concatenate(([0], np.cumsum(lengths)))
    with OfficialImzMLWriter(
        str(path), mode=mode, mz_dtype=np.float64, intensity_dtype=np.float32
    ) as writer:
        for index, coordinate in enumerate(coordinates):
            start, end = offsets[index : index + 2]
            writer.addSpectrum(
                mz if mode == "continuous" else mz[start:end],
                intensity[start:end],
                tuple(coordinate),
            )
    return path, mz, intensity, lengths


@pytest.mark.parametrize("mode", ["continuous", "processed"])
@pytest.mark.parametrize("masked", [False, True])
def test_imzml_export_and_whole_file_coordinate_base(tmp_path, mode, masked):
    coordinates = np.array([[0, 0, 1], [1, 1, 1], [2, 1, 1], [1, 2, 1], [2, 2, 1]])
    source_path, mz, intensity, lengths = _write_imzml_source(
        tmp_path, mode, coordinates
    )
    with MSDataManagerImzML(
        filepath=str(source_path), target_locs=([1, 1], [2, 2])
    ) as source:
        source.preload_hook()
        mask = None
        indices = np.array([1, 2, 3, 4])
        if masked:
            # Source writer reports maxima for a zero-based file; set its canvas explicitly.
            source.ms.meta.max_count_of_pixels_x = 3
            source.ms.meta.max_count_of_pixels_y = 3
            mask = np.zeros((3, 3), dtype=bool)
            mask[1, 2] = mask[2, 2] = True
            indices = np.array([2, 4])
        path = export_imzml(source, tmp_path / "out.imzML", batch_size=1, mask=mask)
        _assert_spectra(path, coordinates, mz, intensity, lengths, indices, mode)


# ---------------------------------------------------------------------------
# Input validation and output cleanup
# ---------------------------------------------------------------------------


def test_ion_empty_mask_precedes_temp_and_writer(tmp_path, monkeypatch):
    with _write_ion_source(tmp_path) as source:

        def unexpected(*args, **kwargs):
            raise AssertionError("writing must not start")

        monkeypatch.setattr(
            "massflow.tools.export_imzml._create_tiled_temp", unexpected
        )
        monkeypatch.setattr(ImzMLWriter, "__init__", unexpected)
        with pytest.raises(ValueError, match="mask does not retain"):
            export_imzml(source, tmp_path / "new" / "out.imzML", mask=np.zeros((5, 6)))
        assert not (tmp_path / "new").exists()


@pytest.mark.parametrize("mode", ["continuous", "processed"])
def test_spectra_empty_selection_cleanup(tmp_path, mode):
    source, _, _, _ = _write_spectra_source(tmp_path, mode)
    with source:
        with pytest.raises(ValueError, match="mask does not retain"):
            export_imzml(
                source, tmp_path / "out.imzML", mask=np.zeros((5, 6)), batch_size=2
            )
        assert not (tmp_path / "out.imzML").exists()
        assert not (tmp_path / "out.ibd").exists()
        source.target_locs = ([99, 99], [100, 100])
        with pytest.raises(ValueError, match="source selection contains no spectra"):
            export_imzml(source, tmp_path / "out.imzML")
        assert not (tmp_path / "out.imzML").exists()


def test_imzml_empty_target_selection_cleanup(tmp_path):
    source_path, _, _, _ = _write_imzml_source(tmp_path, "continuous")
    with MSDataManagerImzML(
        filepath=str(source_path), target_locs=([99, 99], [100, 100])
    ) as source:
        source.preload_hook()
        with pytest.raises(ValueError, match="source selection contains no spectra"):
            export_imzml(source, tmp_path / "out.imzML")
        assert not (tmp_path / "out.imzML").exists()
        assert not (tmp_path / "out.ibd").exists()


def test_preflight_preserves_existing_outputs(tmp_path):
    with _write_ion_source(tmp_path) as source:
        path = tmp_path / "out.imzML"
        path.write_bytes(b"existing XML")
        path.with_suffix(".ibd").write_bytes(b"existing binary")
        with pytest.raises(ValueError, match="mask must have shape"):
            export_imzml(source, path, mask=np.ones((1, 1)))
        assert path.read_bytes() == b"existing XML"
        assert path.with_suffix(".ibd").read_bytes() == b"existing binary"
        with pytest.raises(ValueError, match="mask must have shape"):
            export_imzml(source, tmp_path / "new" / "out.imzML", mask=np.ones((1, 1)))
        assert not (tmp_path / "new").exists()


@pytest.mark.parametrize("state", ["unloaded", "ambiguous", "missing_axis"])
def test_imzml_requires_initialized_state(tmp_path, state):
    source_path, _, _, _ = _write_imzml_source(tmp_path, "continuous")
    with MSDataManagerImzML(filepath=str(source_path)) as source:
        if state != "unloaded":
            source.preload_hook()
            if state == "ambiguous":
                source.ms.meta.processed = True
            else:
                source.ms.shared_mz_list = None
        with pytest.raises(
            ValueError,
            match="loaded continuous/processed mode|initialized shared m/z axis",
        ):
            export_imzml(source, tmp_path / "new" / "out.imzML")
        assert not (tmp_path / "new").exists()


@pytest.mark.parametrize("output_suffix", [".imzML", ".ibd"])
def test_imzml_source_cannot_be_overwritten(tmp_path, output_suffix):
    path, _, _, _ = _write_imzml_source(tmp_path, "processed")
    original_xml = path.read_bytes()
    original_binary = path.with_suffix(".ibd").read_bytes()
    with MSDataManagerImzML(filepath=str(path)) as source:
        source.preload_hook()
        with pytest.raises(ValueError, match="must not overwrite"):
            export_imzml(source, path.with_suffix(output_suffix))
    assert path.read_bytes() == original_xml
    assert path.with_suffix(".ibd").read_bytes() == original_binary


@pytest.mark.parametrize("mode", ["continuous", "processed"])
@pytest.mark.parametrize(
    "invalid", ["intensity", "mz", "lengths", "negative_lengths", "coordinates"]
)
def test_invalid_original_batch_is_checked_before_mask(
    tmp_path, monkeypatch, mode, invalid
):
    source, mz, intensity, lengths = _write_spectra_source(tmp_path, mode)
    coordinates = _COORDINATES
    if invalid == "intensity":
        intensity = intensity[:-1]
    elif invalid == "mz":
        mz = mz[:-1]
    elif invalid == "lengths":
        lengths = lengths.reshape(1, -1)
    elif invalid == "negative_lengths":
        lengths = lengths.copy()
        lengths[0] = -1
    else:
        coordinates = coordinates[:-1]
    with source:
        monkeypatch.setattr(
            source,
            "flat_generator",
            lambda **kwargs: iter([(mz, intensity, lengths, coordinates)]),
        )
        with pytest.raises(ValueError, match="intensity_flat|m/z|lengths|coordinates"):
            export_imzml(source, tmp_path / "out.imzML", mask=np.zeros((5, 6)))
        assert not (tmp_path / "out.imzML").exists()
        assert not (tmp_path / "out.ibd").exists()


def test_unfinalized_source_is_rejected(tmp_path):
    with IonImageDataManagerZarr(
        tmp_path / "source.zarr", mode="w", width=2, height=2, temp_dir=tmp_path
    ) as source:
        with pytest.raises(ValueError, match="must be finalized"):
            export_imzml(source, tmp_path / "out.imzML")
        assert not (tmp_path / "out.imzML").exists()


# ---------------------------------------------------------------------------
# Supported source metadata and writer compatibility
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("masked", [False, True])
def test_supported_zarr_metadata_survives_independent_reopen(tmp_path, masked):
    source, _, _, _ = _write_spectra_source(tmp_path, "continuous")
    try:
        source.add_metadata(
            {
                "centroid_spectrum": False,
                "profile_spectrum": True,
                "polarity": "negative",
                "scan_direction": "bottom_up",
                "line_scan_direction": "line_right_left",
                "scan_pattern": "meandering",
                "scan_type": "vertical_line",
            }
        )
    finally:
        source.close()
    with MSDataManagerZarr(filepath=tmp_path / "spectra.zarr") as source:
        mask = None
        if masked:
            mask = np.zeros((5, 6), dtype=bool)
            mask[1, 2] = True
        path = export_imzml(source, tmp_path / "out.imzML", mask=mask)
    with ImzMLParser(str(path)) as parser:
        assert parser.spectrum_mode == "profile"
        assert parser.polarity == "negative"
        assert parser.root is not None
        for accession in ("IMS:1000400", "IMS:1000490", "IMS:1000410", "IMS:1000481"):
            assert (
                parser.root.find(
                    ".//" + _NS + "cvParam[@accession='" + accession + "']"
                )
                is not None
            )
        assert parser.imzmldict["max count of pixels x"] == (3 if masked else 4)
        assert parser.imzmldict["max count of pixels y"] == (2 if masked else 3)
    with MSDataManagerImzML(filepath=str(path)) as reopened:
        reopened.load_head_data()
        assert reopened.ms.meta.profile_spectrum is True
        assert reopened.ms.meta.mask is not None
        assert reopened.ms.meta.mask.shape == ((2, 3) if masked else (3, 4))
        assert reopened.ms.meta.mask.sum() == (1 if masked else 5)


@pytest.mark.parametrize("mode", ["continuous", "processed"])
def test_supported_imzml_metadata_survives_independent_reopen(tmp_path, mode):
    path = tmp_path / "profile.imzML"
    with OfficialImzMLWriter(
        str(path), mode=mode, spec_type="profile", polarity="positive"
    ) as writer:
        for coordinates, intensity in zip(_COORDINATES, _VALUES):
            writer.addSpectrum(_MZ, intensity, tuple(coordinates))
    with MSDataManagerImzML(filepath=str(path)) as source:
        source.preload_hook()
        output_path = export_imzml(source, tmp_path / "out.imzML")
    with ImzMLParser(str(output_path)) as parser:
        assert parser.spectrum_mode == "profile"
        assert parser.polarity == "positive"
        np.testing.assert_array_equal(parser.coordinates, _COORDINATES)
        for index, intensity in enumerate(_VALUES):
            actual_mz, actual_intensity = parser.getspectrum(index)
            np.testing.assert_array_equal(actual_mz, _MZ)
            np.testing.assert_array_equal(actual_intensity, intensity)


def test_batch_writer_without_metadata_keeps_defaults(tmp_path):
    path = tmp_path / "default.imzML"
    with ImzMLWriter(str(path), mode="continuous") as writer:
        writer.swap_flat_out2disk(_MZ, _VALUES.reshape(-1), np.full(5, 3), _COORDINATES)
    _assert_spectra(
        path,
        _COORDINATES,
        _MZ,
        _VALUES.reshape(-1),
        np.full(5, 3),
        np.arange(5),
        "continuous",
    )
    with ImzMLParser(str(path)) as parser:
        assert parser.spectrum_mode == "centroid"
        assert parser.imzmldict["max count of pixels x"] == 4
        assert parser.imzmldict["max count of pixels y"] == 3
        assert parser.metadata is not None
        assert len(parser.metadata.referenceable_param_groups) == 4


@pytest.mark.parametrize("mode", ["continuous", "processed"])
def test_imzml_one_based_mask_and_empty_mask(tmp_path, mode):
    source_path, mz, intensity, lengths = _write_imzml_source(tmp_path, mode)
    with MSDataManagerImzML(filepath=str(source_path)) as source:
        source.preload_hook()
        mask = np.zeros((3, 4), dtype=bool)
        mask[1, 1] = mask[2, 2] = True
        path = export_imzml(source, tmp_path / "out.imzML", mask=mask, batch_size=2)
        _assert_spectra(
            path, _COORDINATES, mz, intensity, lengths, np.array([0, 4]), mode
        )
        with pytest.raises(ValueError, match="mask does not retain"):
            export_imzml(source, tmp_path / "empty.imzML", mask=np.zeros((3, 4)))
        assert not (tmp_path / "empty.imzML").exists()
        assert not (tmp_path / "empty.ibd").exists()


def test_writer_failure_removes_output_and_transpose(tmp_path, monkeypatch):
    with _write_ion_source(tmp_path) as source:

        def fail_write(*args, **kwargs):
            raise OSError("simulated output failure")

        monkeypatch.setattr(ImzMLWriter, "swap_flat_out2disk", fail_write)
        temp_dir = tmp_path / "transpose"
        with pytest.raises(OSError, match="simulated output failure"):
            export_imzml(source, tmp_path / "out.imzML", temp_dir=temp_dir)
        assert not (tmp_path / "out.imzML").exists()
        assert not (tmp_path / "out.ibd").exists()
        assert not list(temp_dir.iterdir())


@pytest.mark.parametrize("suffix", [".imzML", ".ibd"])
def test_imzml_hardlink_overlap_preserves_source(tmp_path, suffix):
    path, _, _, _ = _write_imzml_source(tmp_path, "processed")
    xml_before = path.read_bytes()
    binary_before = path.with_suffix(".ibd").read_bytes()
    alias = tmp_path / ("alias" + suffix)
    alias.hardlink_to(path.with_suffix(suffix))
    with MSDataManagerImzML(filepath=str(path)) as source:
        source.preload_hook()
        with pytest.raises(ValueError, match="must not overwrite"):
            export_imzml(source, alias)
    assert path.read_bytes() == xml_before
    assert path.with_suffix(".ibd").read_bytes() == binary_before
