
from massflow.module.ms_data_manager_imzml import MSDataManagerImzML
from massflow.module.ms_module import MS
def data_reading_benchmark(mass_data, data_file_path="data/example.imzML"):

    with MSDataManagerImzML(mass_data, filepath=data_file_path) as ms_md:
        ms_md.load_full_data_from_file()
        ms_md.inspect_data()
        for spectrum_data in mass_data:
            _ = spectrum_data.mz_list

def main():
    print("Hello from massflow!")
    mass_data = MS()
    data_reading_benchmark(mass_data)


if __name__ == "__main__":
    main()
