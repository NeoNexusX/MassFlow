"""Tests for the unified MSI metadata model (MetaDataBase / MetaDataImzMl / MetaDataZarr).

Covers:
  - extraction from real imzML fixtures (scope-aware META_INDEX)
  - the historical bugs this refactor fixes (mis-mapped scan_* accessions,
    silent field loss on copy() into a bare MetaDataBase, unreadable
    'flyback' scan pattern)
  - zarr (de)serialization round-trip of the fields formerly duplicated
    across subclasses
  - the ImzMLWriter round-trip for every scan_pattern short key produced
    by populate_from_parser
"""

import numpy as np
import pytest
import zarr
from pyimzml.ImzMLParser import ImzMLParser

from massflow.module.ms_meta_data import (
    META_INDEX,
    OTHER_SECTION,
    MetaDataBase,
    MetaField,
)
from massflow.module.ms_meta_data_imzml import MetaDataImzMl
from massflow.module.ms_meta_data_zarr import MetaDataZarr
from massflow.tools.stream_imzml_writer import ImzMLWriter

EF2193 = "data/Fruitfly_Spleen_IR_MALDESI_20_Negative_ef2193.imzML"
ALIGNED = "data/Fruitfly_Spleen_IR_MALDESI_20_Negative_aligned.imzML"


def _populate(filepath: str) -> MetaDataImzMl:
    parser = ImzMLParser(filepath)
    meta = MetaDataImzMl()
    meta.populate_from_parser(parser)
    return meta


class TestMetaIndexSchema:
    """Guard against META_INDEX drifting away from the declared field schema."""

    def test_every_leaf_target_is_a_declared_metafield(self):
        for accession_map in META_INDEX.values():
            for target in accession_map.values():
                field = target[0] if isinstance(target, tuple) else target
                assert isinstance(
                    getattr(type(MetaDataBase()), field, None), MetaField
                ), f"META_INDEX target {field!r} has no MetaField on MetaDataBase"

    def test_area_keys_have_pyimzml_translation(self):
        from massflow.module.ms_meta_data import AREA_TO_PYIMZML

        assert set(META_INDEX.keys()) == set(AREA_TO_PYIMZML.keys())


class TestPopulateFromParserRealFixtures:
    """Extraction correctness against two real, differently-encoded imzML files."""

    def test_processed_file_file_content_flags(self):
        meta = _populate(EF2193)
        assert meta.processed is True
        assert meta.continuous is None  # mutually exclusive with processed

    def test_processed_file_spectrum_representation(self):
        meta = _populate(EF2193)
        assert meta.profile_spectrum is True
        assert meta.centroid_spectrum is None
        assert meta.ms1_spectrum is True
        assert meta.msn_spectrum is None  # intentionally not auto-extracted

    def test_processed_file_geometry(self):
        meta = _populate(EF2193)
        assert meta.max_count_of_pixels_x == 238
        assert meta.max_count_of_pixels_y == 372
        assert meta.pixel_size_x == pytest.approx(20.0)
        assert meta.pixel_size_y == pytest.approx(20.0)

    def test_processed_file_scan_settings_child_terms(self):
        # This file encodes scan pattern as IMS:1000413 "flyback", a real child
        # term missing from the pre-refactor META_INDEX.
        meta = _populate(EF2193)
        assert meta.scan_pattern == "flyback"
        assert meta.scan_type == "horizontal_line"
        assert meta.scan_direction == "top_down"
        assert meta.line_scan_direction == "line_left_right"

    def test_continuous_file_file_content_flags(self):
        meta = _populate(ALIGNED)
        assert meta.continuous is True
        assert meta.processed is None

    def test_continuous_file_geometry_and_scan_pattern_differ(self):
        meta = _populate(ALIGNED)
        assert meta.max_count_of_pixels_x == 238
        assert meta.max_count_of_pixels_y == 371  # differs from ef2193 (372)
        assert meta.scan_pattern == "one_way"  # differs from ef2193 (flyback)
        # This fixture has no pixel-size cvParam in scanSettingsList.
        assert meta.pixel_size_x is None
        assert meta.pixel_size_y is None

    def test_instrument_model_and_contact_are_not_auto_extracted(self):
        # Category term (MS:1000031) whose value is an unenumerable child
        # accession, and a field intentionally outside META_INDEX.
        for meta in (_populate(EF2193), _populate(ALIGNED)):
            assert meta.instrument_model is None
            assert meta.contact_name == "MassFlow"  # constructor default, untouched


class TestInstrumentComponents:
    """source/analyzer/detector: extracted from componentList's cvParam *names*,
    since (unlike scan_pattern's IMS terms) the PSI-MS source/analyzer/detector
    vocabulary is too open-ended to enumerate as accession->canonical-key pairs.
    Built from a minimal hand-written mzML snippet rather than a real fixture,
    since instrumentConfigurationList only needs source/analyzer/detector's
    componentList children, not a full spectrum-bearing file.
    """

    INSTRUMENT_CONFIG_XML = """<mzML xmlns="http://psi.hupo.org/ms/mzml">
      <fileDescription>
        <fileContent/>
      </fileDescription>
      <instrumentConfigurationList count="1">
        <instrumentConfiguration id="IC1">
          <cvParam accession="MS:1001534" cvRef="MS" name="Bruker Daltonics flex series"/>
          <componentList count="3">
            <source order="1">
              <cvParam accession="MS:1000075" cvRef="MS" name="matrix-assisted laser desorption ionization"/>
            </source>
            <analyzer order="2">
              <cvParam accession="MS:1000084" cvRef="MS" name="time-of-flight"/>
            </analyzer>
            <detector order="3">
              <cvParam accession="MS:1000114" cvRef="MS" name="microchannel plate detector"/>
            </detector>
          </componentList>
        </instrumentConfiguration>
      </instrumentConfigurationList>
    </mzML>"""

    def _populate_from_snippet(self) -> MetaDataImzMl:
        import xml.etree.ElementTree as ET
        from types import SimpleNamespace

        from pyimzml.metadata import Metadata

        root = ET.fromstring(self.INSTRUMENT_CONFIG_XML)
        fake_parser = SimpleNamespace(metadata=Metadata(root))
        meta = MetaDataImzMl()
        meta.populate_from_parser(fake_parser)
        return meta

    def test_source_analyzer_detector_are_extracted_by_name(self):
        meta = self._populate_from_snippet()
        assert meta.source == "matrix-assisted laser desorption ionization"
        assert meta.analyzer == "time-of-flight"
        assert meta.detector == "microchannel plate detector"

    def test_instrument_model_is_still_not_auto_extracted(self):
        # The direct cvParam on <instrumentConfiguration> (the model) is a
        # separate, still-unhandled case from componentList's children.
        meta = self._populate_from_snippet()
        assert meta.instrument_model is None

    def test_instrument_configuration_list_is_a_section_parallel_to_scan_settings(self):
        meta = self._populate_from_snippet()
        assert "instrumentConfigurationList" in meta.sections()
        assert meta.get_section("instrumentConfigurationList") == {
            "source": "matrix-assisted laser desorption ionization",
            "analyzer": "time-of-flight",
            "detector": "microchannel plate detector",
        }

    def test_no_instrument_configuration_is_a_harmless_noop(self):
        from types import SimpleNamespace

        empty_metadata = SimpleNamespace(
            file_description=None,
            referenceable_param_groups={},
            scan_settings={},
            instrument_configurations={},
        )
        meta = MetaDataImzMl()
        meta.populate_from_parser(SimpleNamespace(metadata=empty_metadata))
        assert meta.source is None
        assert meta.get_section("instrumentConfigurationList") == {}


class TestGetSection:
    """Section-based access: META_INDEX segments plus the virtual OTHER_SECTION."""

    def test_sections_lists_index_keys_plus_other(self):
        meta = MetaDataBase()
        assert meta.sections() == [*META_INDEX, OTHER_SECTION]

    def test_populated_fields_grouped_by_section(self):
        meta = _populate(EF2193)
        assert meta.get_section("fileDescription") == {"processed": True}
        assert meta.get_section("referenceableParamGroupList") == {
            "profile_spectrum": True,
            "ms1_spectrum": True,
        }
        assert meta.get_section("scanSettingsList") == {
            "max_count_of_pixels_x": 238,
            "max_count_of_pixels_y": 372,
            "pixel_size_x": pytest.approx(20.0),
            "pixel_size_y": pytest.approx(20.0),
            "scan_pattern": "flyback",
            "scan_type": "horizontal_line",
            "scan_direction": "top_down",
            "line_scan_direction": "line_left_right",
        }

    def test_category_field_appears_once_despite_many_accessions(self):
        # scan_direction is the target of five accessions; it must collapse
        # to a single entry holding the canonical short key.
        meta = MetaDataBase()
        meta.scan_direction = "no_direction"
        section = meta.get_section("scanSettingsList")
        assert section == {"scan_direction": "no_direction"}

    def test_unpopulated_and_none_fields_are_omitted(self):
        meta = MetaDataBase()
        assert meta.get_section("fileDescription") == {}
        meta.continuous = None  # allow_none field explicitly set to None
        assert meta.get_section("fileDescription") == {}

    def test_other_section_holds_fields_outside_meta_index(self):
        meta = MetaDataBase(name="sect_test", instrument_model="Orbitrap")
        other = meta.get_section(OTHER_SECTION)
        assert other["name"] == "sect_test"
        assert other["instrument_model"] == "Orbitrap"
        assert other["contact_name"] == "MassFlow"
        assert other["version"] == 1.0
        # Indexed fields must never leak into "other".
        meta.scan_pattern = "flyback"
        assert "scan_pattern" not in meta.get_section(OTHER_SECTION)

    def test_unknown_section_raises_key_error_listing_valid_names(self):
        meta = MetaDataBase()
        with pytest.raises(KeyError, match="scanSettingsList"):
            meta.get_section("nonexistent")


class TestCopyNoLongerDropsFields:
    """Regression test for the pre-refactor bug: copying into a bare MetaDataBase
    silently discarded fields declared only on the imzML/zarr subclasses."""

    def test_copy_into_bare_base_preserves_previously_subclass_only_fields(self):
        src = MetaDataImzMl()
        src.ms1_spectrum = True
        src.instrument_model = "LTQ FT Ultra"
        src.scan_pattern = "flyback"
        src.absolute_position_offset_x = 12.5

        dst = MetaDataBase()  # mirrors mass_spectrum_set.py's default meta
        dst.copy(src)

        assert dst.ms1_spectrum is True
        assert dst.instrument_model == "LTQ FT Ultra"
        assert dst.scan_pattern == "flyback"
        assert dst.absolute_position_offset_x == pytest.approx(12.5)

    def test_copy_respects_exclude_keys_by_default(self):
        src = MetaDataImzMl()
        src.continuous = True
        src.spectrum_count_num = 999
        dst = MetaDataBase()
        dst.copy(src)
        assert dst.continuous is None
        assert dst.spectrum_count_num is None

    def test_copy_exclude_false_includes_excluded_keys(self):
        src = MetaDataImzMl()
        src.continuous = True
        src.spectrum_count_num = 999
        dst = MetaDataBase()
        dst.copy(src, exclude=False)
        assert dst.continuous is True
        assert dst.spectrum_count_num == 999

    def test_copy_none_source_warns_and_noops(self):
        dst = MetaDataBase()
        dst.name = "untouched"
        dst.copy(None)
        assert dst.name == "untouched"


class TestZarrRoundTrip:
    """Fields uplifted from MetaDataZarr to MetaDataBase must still (de)serialize."""

    def test_roundtrip_preserves_uplifted_fields(self):
        meta = MetaDataZarr(
            name="rt_test",
            max_count_of_pixels_x=10,
            max_count_of_pixels_y=20,
            pixel_size_x=5.0,
            pixel_size_y=5.0,
        )
        meta.ms1_spectrum = True
        meta.instrument_model = "Orbitrap"
        meta.scan_pattern = "flyback"
        meta.continuous = True
        meta.contact_name = "Test Lab"

        root = zarr.group()
        meta.to_zarr(root)
        restored = MetaDataZarr.from_zarr(root)

        assert restored.name == "rt_test"
        assert restored.max_count_of_pixels_x == 10
        assert restored.pixel_size_x == pytest.approx(5.0)
        assert restored.ms1_spectrum is True
        assert restored.instrument_model == "Orbitrap"
        assert restored.scan_pattern == "flyback"
        assert restored.continuous is True
        assert restored.contact_name == "Test Lab"

    def test_roundtrip_mask(self):
        meta = MetaDataZarr(name="mask_test")
        meta.mask = np.zeros((3, 4), dtype=bool)
        root = zarr.group()
        meta.to_zarr(root)
        restored = MetaDataZarr.from_zarr(root)
        assert restored.mask.shape == (3, 4)


class TestWriterScanPatternRoundTrip:
    """Every short key populate_from_parser can produce must be writable by
    ImzMLWriter's obo_codes table (the flyback/no_direction fix)."""

    @pytest.mark.parametrize(
        "pattern",
        ["meandering", "one_way", "random_access", "flyback"],
    )
    def test_scan_pattern_short_key_is_writable(self, tmp_path, pattern):
        out = tmp_path / f"writer_{pattern}"
        writer = ImzMLWriter(str(out), mode="processed", scan_pattern=pattern)
        writer.add_spectrum(np.array([100.0, 200.0]), np.array([10.0, 20.0]), (1, 1))
        writer.close()  # would raise KeyError in obo_codes lookup if unmapped
        assert (tmp_path / f"writer_{pattern}.imzML").exists()

    @pytest.mark.parametrize(
        "direction",
        ["bottom_up", "top_down", "left_right", "right_left", "no_direction"],
    )
    def test_scan_direction_short_key_is_writable(self, tmp_path, direction):
        out = tmp_path / f"writer_dir_{direction}"
        writer = ImzMLWriter(str(out), mode="processed", scan_direction=direction)
        writer.add_spectrum(np.array([100.0]), np.array([10.0]), (1, 1))
        writer.close()
        assert (tmp_path / f"writer_dir_{direction}.imzML").exists()


class TestFullPipelineRealFileLoadThenWrite:
    """End-to-end: real file -> MSDataManagerImzML -> populate_from_parser ->
    the extracted scan_pattern feeds ImzMLWriter without crashing."""

    def test_ef2193_flyback_value_flows_through_writer(self, tmp_path):
        from massflow.data_manager import MSDataManagerImzML

        dm = MSDataManagerImzML(filepath=EF2193)
        dm.load_head_data()
        assert dm.ms.meta.scan_pattern == "flyback"

        out = tmp_path / "roundtrip_from_real_meta"
        writer = ImzMLWriter(
            str(out), mode="processed", scan_pattern=dm.ms.meta.scan_pattern
        )
        writer.add_spectrum(np.array([100.0, 150.0]), np.array([1.0, 2.0]), (1, 1))
        writer.close()
        assert (tmp_path / "roundtrip_from_real_meta.imzML").exists()
