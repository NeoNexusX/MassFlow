
from massflow.module.ms_data_manager_imzml import MSDataManagerImzML
from massflow.module.ms_module import MS
from massflow.logger import get_logger

def main():
    print("Hello from massflow!")
    FILE_PATH = "data/example.imzML"
    ms = MS()
    with MSDataManagerImzML(ms=ms, target_locs=[(1, 1), (50, 50)], filepath=FILE_PATH) as manager:
        manager.load_full_data_from_file()
        manager.inspect_data()
        ms.plot_ms_mask()

if __name__ == "__main__":
    main()