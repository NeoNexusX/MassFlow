"""
ImzML-specific metadata model.

Author: MassFlow Development Team Bionet/dre,NeoNeuxs
License: See LICENSE file in project root
"""
from massflow.module.ms_meta_data import MetaDataBase, MetaField
from massflow.tools.logger import get_logger
logger = get_logger("module")

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
