# module/msi_data_manager_imzml.py

import numpy as np
from pyimzml.ImzMLParser import ImzMLParser, getionimage
import sys

from .msi_data_manager import MSIDataManager
from .msi_module import MSIBaseModule


class MSIDataManagerImzML(MSIDataManager):
    """
    MSI Data Manager for .imzML files.

    Handles reading .imzML files, filtering by m/z range,
    and loading data into the MSI domain model.
    """

    def __init__(self,
                 msi,
                 target_mz_range=None,
                 filepath=None,
                 mz_tolerance=0.1):
        """
        Initialize the ImzML data manager.

        Args:
            msi: MSI object to populate
            target_mz_range: tuple of (min_mz, max_mz) or None for all
            filepath: path to .imzML file
            mz_tolerance: tolerance for ion image extraction (default: 0.1)
        """
        super().__init__(msi, target_mz_range, filepath)
        self.current_num = 0
        self.mz_tolerance = mz_tolerance

    def load_full_data_from_file(self):
        """
        Override the parent method to load data specifically from .imzML.

        If the file is not .imzML, it falls back to the parent's logic
        (which handles .h5 and .msi).
        """
        if self.filepath and self.filepath.lower().endswith('.imzml'):
            self.load_data_from_imzml()
        else:
            print(f"Filepath {self.filepath} is not .imzML. Attempting to load with parent class...", file=sys.stderr)
            super().load_full_data_from_file()


    def load_data_from_imzml(self):
        """
        Main logic for parsing and loading data from an .imzML file.
        """
        print(f"Parsing imzML file: {self.filepath}")

        # 1. Parse the .imzML file
        p = ImzMLParser(self.filepath)

        # 2. Get image dimensions from imzmldict
        max_x = p.imzmldict.get("max count of pixels x", 0)
        max_y = p.imzmldict.get("max count of pixels y", 0)

        if max_x == 0 or max_y == 0:
            print(f"Error: Could not determine image dimensions from imzML file.", file=sys.stderr)
            print(f"imzmldict: {p.imzmldict}", file=sys.stderr)
            return

        shape = (max_y, max_x)
        print(f"Detected image shape: {shape}")

        # Create a mask (imzML doesn't always have a pre-defined mask)
        mask = np.ones(shape, dtype=np.uint8)
        self._msi.meta_mask = mask

        # 3. Collect all unique m/z values from all spectra
        print("Collecting m/z values from all spectra...")
        all_mzs = set()

        # Sample a few spectra to get representative m/z values
        # For large datasets, sampling is more efficient
        sample_size = min(100, len(p.coordinates))
        sample_indices = np.linspace(0, len(p.coordinates) - 1, sample_size, dtype=int)

        for idx in sample_indices:
            mzs, _ = p.getspectrum(idx)
            all_mzs.update(mzs)

        # Convert to sorted array
        mzs_all = np.array(sorted(all_mzs))
        print(f"Found {len(mzs_all)} unique m/z values from sampled spectra")

        # 4. Filter m/z values by target range
        if self.target_mz_range is not None:
            low, high = self.target_mz_range
            mz_mask = (mzs_all >= low) & (mzs_all <= high)
            mzs_to_load = mzs_all[mz_mask]
            print(f"Filtered to {len(mzs_to_load)} m/z values within range {self.target_mz_range}")
        else:
            mzs_to_load = mzs_all
            print(f"Loading all {len(mzs_to_load)} m/z values (this may take a while!)")

        # 5. Set metadata and allocate data matrix
        self._msi.meta_mz_num = len(mzs_to_load)

        if self._msi.meta_mz_num == 0:
            print("No m/z values to load. Aborting.")
            return

        # Allocate the 3D numpy array (mz_num, height, width)
        self._msi.allocate_data_from_meta(dtype=np.float32)

        # 6. Load ion images for each m/z value
        print(f"Loading {len(mzs_to_load)} ion images...")

        for i, mz_value in enumerate(mzs_to_load):
            # Use getionimage to extract the 2D image for this m/z
            msi_image = getionimage(p, mz_value, tol=self.mz_tolerance)

            # Generate base_mask if needed (non-zero pixels)
            base_mask = None
            if self._msi.meta_need_base_mask:
                base_mask = np.where(msi_image > 0, 1, 0).astype(np.uint8)

            # Fill the 3D data matrix
            self._msi.data[self.current_num, :, :] = msi_image

            # Add a slice to the MSI queue
            self._msi.add_msi_slice(
                MSIBaseModule(
                    mz=mz_value,
                    msroi=self._msi.data[self.current_num],  # Reference the slice
                    base_mask=base_mask
                )
            )

            if (i + 1) % 100 == 0 or i == 0:  # Print progress
                print(f'Loaded slice {i + 1}/{self._msi.meta_mz_num} (m/z {mz_value:.4f})')

            self.current_num += 1

        print(f"Finished loading {self.filepath}")
        print(f"Total slices loaded: {self.current_num}")