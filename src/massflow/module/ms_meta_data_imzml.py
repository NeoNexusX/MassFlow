"""
ImzML-specific metadata model.

Author: MassFlow Development Team Bionet/dre,NeoNeuxs
License: See LICENSE file in project root
"""
from typing import Any, Dict, Iterator, Tuple

from pyimzml.metadata import ParamGroup

from massflow.module.ms_meta_data import MetaDataBase, MetaField
from massflow.tools.logger import get_logger
logger = get_logger("massflow.module.meta_data_imzml")

class MetaDataImzMl(MetaDataBase):
    """ImzML metadata wrapper that loads and caches frequently used fields."""

    spectrum_count_num = MetaField()
    absolute_position_offset_x = MetaField()
    absolute_position_offset_y = MetaField()
    instrument_model = MetaField()
    ms1_spectrum = MetaField()
    msn_spectrum = MetaField()

    def __init__(
        self,
        name="imzML_data",
        version=1.0,
        absolute_position_offset_x=None,
        absolute_position_offset_y=None,
        centroid_spectrum=None,
        profile_spectrum=None,
        ms1_spectrum=None,
        msn_spectrum=None,
        instrument_model=None,
        spectrum_count_num=None,
        mask=None,
        pixel_size_x=None,
        pixel_size_y=None,
        max_count_of_pixels_x=None,
        max_count_of_pixels_y=None,
    ):
        """Initialize the metadata object with file path."""
        super().__init__(
            name=name,
            version=version,
            mask=mask,
            pixel_size_x=pixel_size_x,
            pixel_size_y=pixel_size_y,
            max_count_of_pixels_x=max_count_of_pixels_x,
            max_count_of_pixels_y=max_count_of_pixels_y,
        )

        self.absolute_position_offset_x = absolute_position_offset_x
        self.absolute_position_offset_y = absolute_position_offset_y
        self.centroid_spectrum = centroid_spectrum
        self.profile_spectrum = profile_spectrum
        self.ms1_spectrum = ms1_spectrum
        self.msn_spectrum = msn_spectrum
        self.instrument_model = instrument_model
        self.spectrum_count_num = spectrum_count_num

    def copy(self, other: MetaDataBase):
        """Copy specified metadata from other MetaDataImzMl object to self."""

        if other is None:
            logger.warning("Source MS has no metadata to copy or meta is none.")
            return

        exclude_keys = {
            "mask",
            "filepath",
            "spectrum_count_num",
            "continuous",
            "processed",
            "centroid_spectrum",
            "profile_spectrum",
        }

        # Iterate through all attributes and copy if they exist in self and are not excluded
        for key, value in other.items():
            if key in exclude_keys:
                continue

            if hasattr(self, key):
                setattr(self, key, value)

        logger.info("Metadata copy completed.")

    def populate_from_parser(self, parser: Any) -> None:
        """Batch-extract metadata from parser and fill fields defined by `meta_index`."""
        if parser is None:
            logger.error("Parser is not initialized. Please set parser or filepath first.")
            raise ValueError("Parser is not initialized. Please set parser or filepath first.")

        values_by_accession = self._collect_accession_values(parser)
        for accession_id, prop_name in self.meta_index.items():
            param_value = values_by_accession.get(accession_id)
            if param_value is not None:
                setattr(self, prop_name, param_value)

    def _collect_accession_values(self, parser: Any) -> Dict[str, Any]:
        """Traverse all metadata areas once and collect values for required accession ids."""
        required_accessions = set(self.meta_index.keys())
        found: Dict[str, Any] = {}

        search_areas = [
            parser.metadata.file_description,
            parser.metadata.scan_settings,
            parser.metadata.instrument_configurations,
            parser.metadata.samples,
            parser.metadata.softwares,
            parser.metadata.data_processings,
            parser.metadata.referenceable_param_groups,
        ]

        for area in search_areas:
            for accession_id, value in self._iter_area_params(area, required_accessions):
                if accession_id not in found:
                    found[accession_id] = value
                    if len(found) == len(required_accessions):
                        return found

        return found

    def _iter_area_params(self, area: Any, required: set) -> Iterator[Tuple[str, Any]]:
        """Yield (accession_id, value) pairs from one metadata area."""
        if area is None:
            return

        if isinstance(area, ParamGroup):
            yield from self._iter_param_group(area, required)
            return

        if isinstance(area, dict):
            for param_group in area.values():
                if isinstance(param_group, ParamGroup):
                    yield from self._iter_param_group(param_group, required)

    @staticmethod
    def _iter_param_group(param_group: ParamGroup, required: set) -> Iterator[Tuple[str, Any]]:
        """Yield all accession/value entries from one ParamGroup."""
        for accession_id in required:
            if accession_id in param_group:
                yield accession_id, param_group[accession_id]
