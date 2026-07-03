"""
ImzML-specific metadata model.

Author: MassFlow Development Team Bionet/dre,NeoNeuxs
License: See LICENSE file in project root
"""
from massflow.module.ms_meta_data import MetaDataBase
from massflow.tools.logger import get_logger
logger = get_logger("massflow.module.meta_data_imzml")

class MetaDataImzMl(MetaDataBase):
    """ImzML metadata wrapper.

    Field schema and `populate_from_parser` extraction live on `MetaDataBase`;
    this subclass only supplies the imzML-specific default name.
    """

    def __init__(self, name="imzML_data", **kwargs):
        super().__init__(name=name, **kwargs)
