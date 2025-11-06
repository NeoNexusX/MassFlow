import sys
import os
from pathlib import Path
import numpy as np

# Allow importing modules from the project root directory
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from module.ms_module import MS
from module.ms_data_manager_imzml import MSDataManagerImzML
from ms_preprocess import MSIPreprocessor


def main():
    # Input file (modify as needed)
    file_path = str(Path("Dataset/neg-gz4.imzML"))

    # Load imzML into MS queue
    ms = MS()
    ms_dm = MSDataManagerImzML(ms, filepath=file_path)
    ms_dm.load_full_data_from_file()

    if len(ms) == 0:
        print("No spectra loaded. Check if the imzML path is correct.")
        return

    # Take the first spectrum for denoising
    spectrum = ms[0]

    # Select simple and robust moving average (only depends on numpy)
    method = "ma"
    params = {"window": 7}  # Adjust window size as needed (will be automatically converted to odd)

    try:
        denoised = MSIPreprocessor.noise_reduction(spectrum, method=method, **params)
    except ImportError as e:
        print(f"Failed to import required libraries: {e}")
        print("If you switch to 'gaussian'/'savgol'/'wavelet', install: scipy, pywavelets")
        return
    except Exception as e:
        print(f"Noise reduction failed: {e}")
        import traceback
        traceback.print_exc()
        return

    # Save denoising results (intensity and corresponding m/z)
    out_base = "neg-gz4_ma_denoised"
    inten_den = np.array(denoised.intensity)
    mz_den = np.array(denoised.mz_list)

    np.save(f"{out_base}_intensity.npy", inten_den)
    np.save(f"{out_base}_mz.npy", mz_den)

    print(f"Saved denoised intensity: {out_base}_intensity.npy")
    print(f"Saved denoised m/z: {out_base}_mz.npy")
    print("Done.")


if __name__ == "__main__":
    main()