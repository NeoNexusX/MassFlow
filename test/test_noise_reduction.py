import pytest
from massflow.module.ms_module import MS
from massflow.module.ms_data_manager_imzml import MSDataManagerImzML
from massflow.logger import get_logger

logger = get_logger("__test_tools__")

@pytest.fixture(scope="class")
def pre_load_data(data_file_path="data/example.imzML"):
    #pre load data
    mass_data = MS()

    with MSDataManagerImzML(mass_data, filepath=data_file_path) as ms_md:
        ms_md.load_full_data_from_file()
        ms_md.inspect_data()

        # Pre-load data
        for spectrum_data in mass_data:
            _ = spectrum_data.mz_list
        logger.info("data pre-load finished!")

    yield mass_data

    # Cleanup code can be added here if needed
    del ms_md
    del mass_data

class TestNoiseReduction:

    def test_noise_reduction_time(self, pre_load_data):
        pass
        # Your test code here

    def test_noise_reduction_snr(self, pre_load_data):
        # Your test code here
        pass