import os
import time
import numpy as np
import pytest
from massflow.logger import get_logger

from massflow.module.mass_spectrum_set import MassSpectrumSet
from massflow.module.ms_data_manager_imzml import MSDataManagerImzML
from massflow.tools.plot import plot_two_together, plot_two_individual, plot_spectrum
from massflow.preprocess.ms_preprocess import MSIPreprocessor

FILE_PATH = 'data/neg-gz4.imzML'
ms = MassSpectrumSet()
ms_md = MSDataManagerImzML(ms, filepath=FILE_PATH)
ms_md.load_full_data_from_file()
sp = ms[0]

logger = get_logger("test_plot")
MZ_RANGE = (500.0, 550.0)
LOG_DIR = "logs"

denoised_spectrum = MSIPreprocessor.noise_reduction_spectrum(
    data=sp,
    method="ma",
    window=7  # Window size for moving average
)

class TestPlot:

    def test_plot_two_together_single_spectrum(self):
        os.makedirs(LOG_DIR, exist_ok=True)
        save_path = os.path.join(LOG_DIR, "plot_together_single_line.png")
        logger.info("Testing plot_two_together with line mode")
        plot_two_together(target=denoised_spectrum, base=sp, save_path=save_path, mz_range=MZ_RANGE, intensity_range=(0.0, 1.5), metrics_box=True, plot_mode=["line","line"])
        assert os.path.exists(save_path)

    def test_plot_two_individual_single_spectrum(self):
        os.makedirs(LOG_DIR, exist_ok=True)
        save_path = os.path.join(LOG_DIR, "plot_two_individual_single_line.png")
        logger.info("Testing plot_two_individual with line mode")
        plot_two_individual(target=denoised_spectrum, base=sp, save_path=save_path, mz_range=MZ_RANGE, intensity_range=(0.0, 1.5), metrics_box=True, plot_mode=["line","line"])
        assert os.path.exists(save_path)

    def test_plot_spectrum_single_spectrum(self):
        os.makedirs(LOG_DIR, exist_ok=True)
        save_path = os.path.join(LOG_DIR, "plot_single_line.png")
        logger.info("Testing plot_spectrum single with line mode")
        plot_spectrum(base=denoised_spectrum, save_path=save_path, mz_range=MZ_RANGE, metrics_box=False, intensity_range=(0.0, 1.5), title_suffix="single", plot_mode=["line"])
        assert os.path.exists(save_path)