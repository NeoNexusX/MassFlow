import numpy as np
from massflow.module.mass_spectrum_set import MassSpectrumSet
from massflow.module.ms_data_manager_imzml import MSDataManagerImzML
import pytest
import time
from functools import partial
from massflow.logger import get_logger

FILE_PATH = "data/neg-gz4.imzML"

def run_crop(mass_data: MassSpectrumSet, mz_range, mode='update', batch_size=256, max_threads=-1):
    with MSDataManagerImzML(mass_data, filepath=FILE_PATH) as manager:
        manager.ms = mass_data
        result = manager.update_target_mz(mz_range=mz_range, batch_size=batch_size, max_threads=max_threads, mode=mode)
        return result


logger = get_logger("test_crop_range")

def pre_load_data_setup(data_file_path=FILE_PATH):
    mass_data = MassSpectrumSet()
    with MSDataManagerImzML(mass_data, filepath=data_file_path) as ms_md:
        ms_md.load_full_data_from_file()
        # Pre-load intensity to simulate in-memory processing benchmark
        for spectrum_data in mass_data:
            _ = spectrum_data.intensity
        logger.info("data pre-load finished!")
    return (mass_data,), {}

class TestCropRangeBenchmark:

    @pytest.mark.parametrize("mode", ["update","new"])
    @pytest.mark.benchmark(timer=time.perf_counter)
    def test_update_target_mz_benchmark(self, benchmark, mode):
        mz_range = (400.0, 500.0)
        logger.info(f"Testing update_target_mz (mode={mode})")
        benchmark.pedantic(partial(run_crop, mz_range=mz_range, mode=mode),
                           rounds=1,
                           setup=pre_load_data_setup,
                           iterations=1,
                           warmup_rounds=0)
        logger.info("Finished update_target_mz benchmark")


    # def test_crop_range_functionality(self):
    #     """
    #     Verify cropping occurs with a valid in-range window for both 'update' and 'new' modes.
    #     """
    #     mz_range = (500.0, 600.0)

    #     args, _ = pre_load_data_setup()
    #     ms_update = args[0]
    #     orig_lens_update = [len(sp.mz_list) for sp in ms_update]

    #     run_crop(ms_update, mz_range, mode='update')

    #     any_reduced = False
    #     for i, sp in enumerate(ms_update):
    #         cur_len = len(sp.mz_list)
    #         if cur_len < orig_lens_update[i]:
    #             any_reduced = True
    #         # Skip empties marked by sentinel [-1]
    #         if isinstance(sp.mz_list, np.ndarray) and sp.mz_list.size == 1 and sp.mz_list[0] == -1:
    #             continue
    #         # Non-empty spectra must lie within the crop window
    #         assert sp.mz_list[0] >= mz_range[0]
    #         assert sp.mz_list[-1] <= mz_range[1]
    #         assert cur_len <= orig_lens_update[i]
    #     assert any_reduced, "At least one spectrum should be reduced after cropping in update mode"

    #     # New mode: ensure returned object has spectra within window and same count
    #     args, _ = pre_load_data_setup()
    #     ms_new_src = args[0]
    #     orig_lens_new = [len(sp.mz_list) for sp in ms_new_src]

    #     new_ms = run_crop(ms_new_src, mz_range, mode='new')
    #     assert new_ms is not ms_new_src, "Should return a new MS object"
    #     assert len(new_ms) == len(ms_new_src), "New MS should have same number of spectra"

    #     for i, sp in enumerate(new_ms):
    #         cur_len = len(sp.mz_list)
    #         # Skip empties marked by sentinel [-1]
    #         if isinstance(sp.mz_list, np.ndarray) and sp.mz_list.size == 1 and sp.mz_list[0] == -1:
    #             continue
    #         assert sp.mz_list[0] >= mz_range[0]
    #         assert sp.mz_list[-1] <= mz_range[1]
    #         assert cur_len <= orig_lens_new[i]

    # def test_invalid_range_sets_intensity_sentinel_update(self):
    #     args, _ = pre_load_data_setup()
    #     ms_update = args[0]
    #     min_mz = float(np.min(ms_update[0].mz_list))
    #     invalid_range = (min_mz - 100.0, min_mz - 50.0)
    #     run_crop(ms_update, invalid_range, mode='update')
    #     for sp in ms_update:
    #         assert isinstance(sp.intensity, np.ndarray)
    #         assert sp.intensity.size == 1 and sp.intensity[0] == -1, "Intensity must be sentinel -1 for invalid range"

    # def test_invalid_range_sets_intensity_sentinel_new(self):
    #     args, _ = pre_load_data_setup()
    #     ms_src = args[0]
    #     min_mz = float(np.min(ms_src[0].mz_list))
    #     invalid_range = (min_mz - 100.0, min_mz - 50.0)
    #     new_ms = run_crop(ms_src, invalid_range, mode='new')
    #     assert new_ms is not ms_src
    #     for sp in new_ms:
    #         assert isinstance(sp.intensity, np.ndarray)
    #         assert sp.intensity.size == 1 and sp.intensity[0] == -1, "Intensity must be sentinel -1 for invalid range"