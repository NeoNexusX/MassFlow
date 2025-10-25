"""
MSI Data Manager for imzML format.

Handles reading .imzML files and loading data into the MSI domain model.
Supports filtering by m/z range and efficient ion image extraction.
"""

import numpy as np
from pyimzml.ImzMLParser import ImzMLParser, getionimage
import sys
import os

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
                 mz_tolerance=0.1,
                 sample_spectra=100):
        """
        Initialize the ImzML data manager.

        Args:
            msi: MSI object to populate
            target_mz_range: tuple of (min_mz, max_mz) or None for all
            filepath: path to .imzML file
            mz_tolerance: tolerance for ion image extraction (default: 0.1)
            sample_spectra: number of spectra to sample for m/z collection (default: 100)
        """
        super().__init__(msi, target_mz_range, filepath)
        self.mz_tolerance = mz_tolerance
        self.sample_spectra = sample_spectra

    def load_full_data_from_file(self):
        """
        Override the parent method to load data specifically from .imzML.

        If the file is not .imzML, it falls back to the parent's logic.
        """
        if not self.filepath:
            print("Error: No filepath provided.", file=sys.stderr)
            return

        if not os.path.exists(self.filepath):
            print(f"Error: File {self.filepath} does not exist.", file=sys.stderr)
            return

        if self.filepath.lower().endswith('.imzml'):
            self.load_data_from_imzml() #一切正常的话，这个方法里主要就执行这一句
        else:
            print(f"Warning: Filepath {self.filepath} is not .imzML. "
                  "Attempting to load with parent class...", file=sys.stderr)
            super().load_full_data_from_file()

    def load_data_from_imzml(self): #整体的工作流程在这
        """
        Main logic for parsing and loading data from an .imzML file.
        """
        print(f"Parsing imzML file: {self.filepath}")

        try:
            # 1. Parse the .imzML file
            parser = ImzMLParser(self.filepath)
        except Exception as e:
            print(f"Error: Failed to parse imzML file: {e}", file=sys.stderr)
            return

        # 2. Get image dimensions
        shape = self._get_image_shape(parser)
        if shape is None:
            return

        print(f"Detected image shape: {shape}")

        # 3. Create and set mask
        mask = self._create_mask(parser, shape)
        self._msi.meta_mask = mask

        # 4. Collect m/z values
        mzs_to_load = self._collect_mz_values(parser) #给mzs_to_load的返回值是一个列表，里面是通道值
        if mzs_to_load is None or len(mzs_to_load) == 0:
            print("No m/z values to load. Aborting.")
            return

        # 5. Set metadata and allocate data matrix
        self._msi.meta_mz_num = len(mzs_to_load)
        self._msi.allocate_data_from_meta(dtype=np.float32)

        # 6. Load ion images
        self._load_ion_images(parser, mzs_to_load, shape)

        print(f"Finished loading {self.filepath}")
        print(f"Total slices loaded: {self.current_image_num}")

        #同步元数据
        self._msi.update_metadata()

    def _get_image_shape(self, parser):
        """
        Extract image dimensions from the imzML parser.

        Args:
            parser: ImzMLParser instance

        Returns:
            tuple: (height, width) or None if dimensions cannot be determined
        """
        max_x = parser.imzmldict.get("max count of pixels x", 0)
        max_y = parser.imzmldict.get("max count of pixels y", 0)

        if max_x == 0 or max_y == 0:
            print("Error: Could not determine image dimensions from imzML file.",
                  file=sys.stderr)
            print(f"Available metadata: {parser.imzmldict}", file=sys.stderr)
            return None

        return (max_y, max_x)


    def _create_mask(self, parser, shape):
        """
        Create a mask based on actual pixel coordinates.

        Args:
            parser: ImzMLParser instance
            shape: tuple of (height, width)

        Returns:
            numpy.ndarray: 2D binary mask
        """
        mask = np.zeros(shape, dtype=np.uint8)

        try:
            # 使用更快的 NumPy 索引
            # 1. 解包所有坐标 (imzML 坐标是 1-based)
            # 确保 parser.coordinates 是一个列表
            coords = list(parser.coordinates)
            x_indices = [c[0] - 1 for c in coords]
            y_indices = [c[1] - 1 for c in coords]

            # 2. 使用 NumPy 一次性赋值
            mask[y_indices, x_indices] = 1

        except Exception:
            # --- 如果失败，回退到安全循环 ---
            print("Warning: NumPy indexing failed, falling back to iterative mask creation.")
            for x, y, _ in parser.coordinates:
                # imzML 坐标是 1-indexed, convert to 0-indexed
                mask[y - 1, x - 1] = 1


        non_zero_pixels = np.sum(mask)
        print(f"Created mask with {non_zero_pixels} non-zero pixels out of "
              f"{shape[0] * shape[1]} total pixels")

        return mask


    def _collect_mz_values(self, parser): #收集离子通道数
        """
        Collect unique m/z values.
        Tries to get a global m/z list first (Processed Data).
        If fails, falls back to sampling (Continuous Data).
        """
        print("Collecting m/z values...")
        mzs_all = None

        # 1. 尝试“已处理型”路径 (快速)
        if hasattr(parser, 'mzs') and len(parser.mzs) > 0:
            print("Found global m/z list in file metadata.")
            mzs_all = np.array(parser.mzs)

        # 2. 如果快速路径失败，回退到“连续型”抽样路径
        else:
            print("No global m/z list found. Sampling spectra...")
            total_spectra = len(parser.coordinates)
            sample_size = min(self.sample_spectra, total_spectra)

            if sample_size == 0: #采样的像素点不能为0
                print("Error: No spectra coordinates found.", file=sys.stderr)
                return None

            sample_indices = np.linspace(0, total_spectra - 1, sample_size, dtype=int)#生成sample_indices

            all_mzs_set = set()  # 使用新名称避免混淆
            for idx in sample_indices:
                try:
                    mzs, _ = parser.getspectrum(idx)
                    all_mzs_set.update(mzs)
                except Exception as e:
                    print(f"Warning: Failed to read spectrum {idx}: {e}", file=sys.stderr)
                    continue

            # 安全检查
            if not all_mzs_set:
                print("Error: No m/z values could be collected from sampling.", file=sys.stderr)
                return None

            mzs_all = np.array(sorted(all_mzs_set))

        # 3. 对来自任一路径的结果进行过滤
        print(f"Found {len(mzs_all)} unique m/z values")

        if self.target_mz_range is not None:
            low, high = self.target_mz_range
            mz_mask = (mzs_all >= low) & (mzs_all <= high)
            mzs_to_load = mzs_all[mz_mask]
            print(f"Filtered to {len(mzs_to_load)} m/z values within range "
                  f"[{low}, {high}]")
        else:
            mzs_to_load = mzs_all
            print(f"Loading all {len(mzs_to_load)} m/z values")

        return mzs_to_load

    def _load_ion_images(self, parser, mzs_to_load, shape):
        """
        Load ion images for each m/z value.

        Args:
            parser: ImzMLParser instance
            mzs_to_load: numpy.ndarray of m/z values
            shape: tuple of (height, width)
        """
        print(f"Loading {len(mzs_to_load)} ion images...")

        for i, mz_value in enumerate(mzs_to_load):
            try:
                # Extract 2D image for this m/z
                msi_image = getionimage(parser, mz_value, tol=self.mz_tolerance)

                # Handle potential shape mismatch
                if msi_image.shape != shape:
                    print(f"Warning: Image shape {msi_image.shape} differs from "
                          f"expected {shape} for m/z {mz_value:.4f}", file=sys.stderr)
                    # Resize or pad if needed
                    msi_image = self._resize_image(msi_image, shape)

                # Generate base_mask if needed
                base_mask = None
                if self._msi.meta_need_base_mask:
                    base_mask = np.where(msi_image > 0, 1, 0).astype(np.uint8)

                # Fill the 3D data matrix
                self._msi.data[self.current_image_num, :, :] = msi_image

                # Add a slice to the MSI queue
                self._msi.add_msi_slice(
                    MSIBaseModule(
                        mz=mz_value,
                        msroi=self._msi.data[self.current_image_num],
                        base_mask=base_mask
                    )
                )

                # Progress reporting
                if (i + 1) % 100 == 0 or i == 0:
                    print(f'Loaded slice {i + 1}/{len(mzs_to_load)} '
                          f'(m/z {mz_value:.4f})')

                self.current_image_num += 1

            except Exception as e:
                print(f"Error loading m/z {mz_value:.4f}: {e}", file=sys.stderr)
                continue

    @staticmethod
    def _resize_image(image, target_shape):
        """
        Resize or pad image to match target shape.

        Args:
            image: numpy.ndarray, input image
            target_shape: tuple of (height, width)

        Returns:
            numpy.ndarray: resized image
        """
        result = np.zeros(target_shape, dtype=image.dtype)
        h, w = min(image.shape[0], target_shape[0]), min(image.shape[1], target_shape[1])
        result[:h, :w] = image[:h, :w]
        return result