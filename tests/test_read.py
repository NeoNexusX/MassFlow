from massflow.logger import get_logger
logger = get_logger("test_read")

def data_reading_benchmark(mass_data, data_file_path="data/example.imzML"):
    with MSDataManagerImzML(mass_data, filepath=data_file_path) as ms_md:
        ms_md.load_full_data_from_file()
        ms_md.inspect_data()