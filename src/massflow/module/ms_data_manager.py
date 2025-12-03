"""
MS Data Management Module

Provides functions for reading/writing MS data, memory statistics, and visualization.
Supports .h5/.msi files and batch import from directories, filters by m/z range,
and generates merged or split outputs.

Author: MassFlow Development Team Bionet/NeoNexus
License: See LICENSE file in project root
"""
from abc import ABC, abstractmethod
from typing import Tuple
import numpy as np
from massflow.logger import get_logger
from massflow.module.ms_module import MS
from concurrent.futures import ThreadPoolExecutor,ProcessPoolExecutor
import math
logger = get_logger("ms_data_manager")

def _read_spectrum_raw(args: Tuple[str, object, int]) -> Tuple[int, np.ndarray, np.ndarray]:
    """
    子进程执行：根据 ibd_path + index 读取一个 spectrum 的 (mz, intensity)。

    返回:
        (index, mz, intensity)
    """
    ibd_path, reader, index = args
    with open(ibd_path, "rb") as f:
        mz, intensity = reader.read_spectrum_from_file(f, index)
    return index, mz, intensity

def _process_indices_in_threads(args: Tuple[str, object, list[int], int]):
    """
    子进程内部使用多线程读取一批索引。

    Args:
        ibd_path (str)
        reader (object)
        indices (list[int])
        thread_workers (int)

    Returns:
        list[tuple[int, np.ndarray, np.ndarray]]
    """
    ibd_path, reader, indices, thread_workers = args

    def _read_one(idx: int):
        return _read_spectrum_raw((ibd_path, reader, idx))

    results: list[tuple[int, np.ndarray, np.ndarray]] = []
    # 线程数不要超过索引数
    n_threads = min(thread_workers, max(1, len(indices)))
    with ThreadPoolExecutor(max_workers=n_threads) as executor:
        for res in executor.map(_read_one, indices):
            results.append(res)
    return results

class MSDataManager(ABC):
    """
    Abstract base class for MS Data Manager.
    """
    def __init__(self,
                 ms: MS,
                 target_mz_range=None,
                 target_locs=None,
                 filepath=None):
        """
        Initialize the MS data manager.

        Args:
            ms (MS): Mass-spectrometry domain model instance.
            target_mz_range (tuple[float, float], optional): (min_mz, max_mz) to filter peaks.
            target_locs (list[tuple], optional): List of (x,y) or (x,y,z) coordinates to load.
            filepath (str, optional): Path to the input file.
        """
        self._ms = ms
        self.target_mz_range = target_mz_range

        # target_locs input verification
        if target_locs is not None:
            # target locs input
            if len(target_locs) <= 1:
                logger.error("target_locs must be non-empty")
                raise ValueError("target_locs must be non-empty")
            elif target_locs[0][0] > target_locs[1][0] or target_locs[0][1] > target_locs[1][1]:
                logger.error("locs must x1<x2,y1<y2")
                raise ValueError("locs must x1<x2,y1<y2")

        # update target_locs
        self.target_locs = target_locs
        self.filepath = filepath
        self.current_spectrum_num = 0

    @property
    def ms(self) -> MS:
        """
        Get the MS object.

        Returns:
            MS: The MS object.
        """
        return self._ms

    @ms.setter
    def ms(self, ms: MS):
        """
        Set the MS object.

        Args:
            ms (MS): Mass-spectrometry domain model instance.
        """
        self._ms = ms

    @abstractmethod
    def load_full_data_from_file(self):
        """
        Load metadata from a file.

        Args:
            filepath (str): Path to the input file.
        """

    def inspect_data(self,inpect_num=2):
        """
        Inspect the data structure of the MSI object.

        Log metadata shapes and queue information, including max/min m/z values,
        queue length, and count of non-empty base masks.
        """
        meta_info = "MS meta data:\r\n"
        meta_info+=(f"  target_mz_range: {self.target_mz_range}\r\n")
        meta_info+=(f"  target_locs: {self.target_locs}\r\n")
        meta_info+=(f"  filepath: {self.filepath}\r\n")
        meta_info+=(f"  current_spectrum_num: {self.current_spectrum_num}\r\n")
        
        # get meta data info
        if self.ms.meta is not None:
            for attr, value in self.ms.meta.items():
                shape = getattr(value, 'shape', None)
                if shape is not None and len(shape) > 0:
                    meta_info+=(f"  meta_{attr}: {shape}\r\n")
                else:
                    meta_info+=(f"  meta_{attr}: {value}\r\n")
            logger.info(meta_info)

        # get base info
        base_info = "MS  information:\r\n"
        pointer4num = 0
        for spectrum in self._ms:
            if pointer4num >= inpect_num:
                break
            base_info+=(f"  MS len: {len(spectrum)}\r\n"
                        f"  MS range: {min(spectrum.mz_list)} - {max(spectrum.mz_list)}\r\n"
                        f"  MS coord: {spectrum.coordinates}\r\n"
                        f"  max and min mz_list: {max(spectrum.mz_list)} - {min(spectrum.mz_list)}\r\n"
                        f"  max intensity: {max(spectrum.intensity)}\r\n\r\n")
            pointer4num += 1
        logger.info(base_info)

    def get_batch_generator(self, batch_size: int = 32, max_workers: int = 32):
        total_length = len(self.ms)

        def _load_slice(slice_data):
            for sp in slice_data:
                _ = sp.intensity

        for i in range(0, total_length, batch_size):
            batch = self.ms[i: i + batch_size]
            batch_len = len(batch)
            if batch_len <= 1 or max_workers <= 1:
                for sp in batch:
                    _ = sp.intensity #type: ignore
                yield batch
                continue

            slice_size = max(1, (batch_len + max_workers - 1) // max_workers)
            with ThreadPoolExecutor(max_workers=min(max_workers, batch_len)) as executor:
                futures = []
                for start in range(0, batch_len, slice_size):
                    end = min(start + slice_size, batch_len)
                    futures.append(executor.submit(_load_slice, batch[start:end]))
                for f in futures:
                    f.result()
            yield batch

    def get_multiprocess_batch_generator(self, batch_size: int = 512, max_workers: int = 4):
        """
        使用多进程按 batch 预加载 spectra 的 mz_list 和 intensity，并按 batch 返回。

        每个 batch 的流程:
            - 根据 batch 的索引范围 [start, end) 构造参数列表；
            - ProcessPoolExecutor 并行调用 _read_spectrum_raw；
            - 主进程根据返回的 index 把 mz/intensity 写回 self.ms[index]；
            - 收集该 batch 对应的 Spectrum 对象列表并 yield。

        Args:
            batch_size (int): 每个 batch 的谱条数。
            max_workers (int): 进程池最大进程数。

        Yields:
            list[SpectrumBaseModule]: 当前 batch 的光谱对象列表（已完成预加载）。
        """
        total = len(self.ms)
        if total == 0:
            logger.info("No spectra to preload in batches.")
            return

        # 假定所有 spectrum 共享同一个 reader / ibd_path（与 preload_all_spectra_multiprocess 相同）
        sample_sp = self.ms[0]
        ibd_path = sample_sp._ibd_path # type: ignore
        reader = sample_sp._reader # type: ignore

        logger.info(
            f"Start multiprocess batch preload: total={total}, "
            f"batch_size={batch_size}, max_workers={max_workers}"
        )

        # 这里选择“每个 batch 建一个进程池”，避免单个长生命周期 executor 里任务过多导致调度不均
        for start in range(0, total, batch_size):
            end = min(start + batch_size, total)
            batch_indices = list(range(start, end))

            # 构造本 batch 的参数
            args_iter = ((ibd_path, reader, idx) for idx in batch_indices)

            # 预加载当前 batch
            with ProcessPoolExecutor(max_workers=max_workers) as executor:
                for idx, mz, inten in executor.map(_read_spectrum_raw, args_iter):
                    sp = self.ms[idx]
                    sp.mz_list = mz.copy() # type: ignore
                    sp.intensity = inten.copy() # type: ignore

            # 预加载完成后，构造并返回本 batch 的 spectrum 列表
            batch = [self.ms[idx] for idx in batch_indices]
            yield batch

        logger.info("Finished multiprocess batch preload for all spectra.")

    def get_multiprocess_batch_generator_hybrid(
            self,
            batch_size: int = 512,
            process_workers: int = 4,
            thread_workers_per_process: int = 64,
    ):
        """
        混合模式：多进程 + 每进程多线程 预加载 spectra 的 mz_list 和 intensity。

        每个 batch 的流程:
            - 将 batch 内的索引平均切成 process_workers 份；
            - 为每一份启动一个进程；
            - 每个进程内部再用 ThreadPoolExecutor 读取该子集索引；
            - 主进程将所有返回结果写回 self.ms，并按顺序 yield batch。

        Args:
            batch_size (int): 每个 batch 的谱条数。
            process_workers (int): 进程池最大进程数。
            thread_workers_per_process (int): 每个进程内部的最大线程数。

        Yields:
            list[SpectrumBaseModule]: 当前 batch 的光谱对象列表（已完成预加载）。
        """
        total = len(self.ms)
        if total == 0:
            logger.info("No spectra to preload in batches (hybrid).")
            return

        sample_sp = self.ms[0]
        ibd_path = sample_sp._ibd_path  # type: ignore
        reader = sample_sp._reader      # type: ignore

        logger.info(
            f"Start HYBRID multiprocess batch preload: total={total}, "
            f"batch_size={batch_size}, process_workers={process_workers}, "
            f"threads_per_proc={thread_workers_per_process}"
        )

        for start in range(0, total, batch_size):
            end = min(start + batch_size, total)
            batch_indices = list(range(start, end))
            batch_len = len(batch_indices)
            if batch_len == 0:
                continue

            # 按进程数切分索引
            n_procs = min(process_workers, batch_len)
            chunk_size = math.ceil(batch_len / n_procs)
            index_chunks: list[list[int]] = [
                batch_indices[i:i + chunk_size]
                for i in range(0, batch_len, chunk_size)
            ]

            # 为每个 chunk 构造参数：在子进程里再用多线程读这一块
            proc_args = [
                (ibd_path, reader, idx_chunk, thread_workers_per_process)
                for idx_chunk in index_chunks
            ]

            results: list[tuple[int, np.ndarray, np.ndarray]] = []
            # 多进程：每个进程内部再用多线程
            with ProcessPoolExecutor(max_workers=n_procs) as proc_pool:
                for chunk_result in proc_pool.map(_process_indices_in_threads, proc_args):
                    # chunk_result 是当前子进程返回的一批 (idx, mz, inten)
                    results.extend(chunk_result)

            # 为了安全起见，对结果按 index 排序再写回（避免乱序）
            results.sort(key=lambda x: x[0])

            for idx, mz, inten in results:
                sp = self.ms[idx]
                sp.mz_list = mz.copy()      # type: ignore
                sp.intensity = inten.copy() # type: ignore

            batch = [self.ms[idx] for idx in batch_indices]
            yield batch

        logger.info("Finished HYBRID multiprocess batch preload for all spectra.")

    def create_ms_meta_mask(self):
        """
        Create a binary occupancy mask for all available spectra coordinates.

        Parameters
        - None

        Returns
        - np.ndarray: A 2D array of shape `(height, width)` where `1` marks a pixel
          with an available spectrum and `0` otherwise.

        Raises
        - ValueError: If required metadata (`max_count_of_pixels_x/y`) is missing.
        """
        if self.ms.meta is None:
            logger.error("MS meta data is None. Please load meta data first.")
            raise ValueError("MS meta data is None. Please load meta data first.")

        if self.ms.meta.max_count_of_pixels_x is None or self.ms.meta.max_count_of_pixels_y is None:
            logger.error("Image dimensions missing in meta data.")
            raise ValueError("Image dimensions missing in meta data.")

        width = int(self.ms.meta.max_count_of_pixels_x)
        height = int(self.ms.meta.max_count_of_pixels_y)

        # Prepare an empty occupancy mask
        mask = np.zeros((height, width), dtype=np.int8)

        # Fill the occupancy mask across all z-planes
        for z_dict in self.ms.coordinate_index.values():
            for x, y_dict in z_dict.items():
                for y in y_dict.keys():
                    ix = int(x)
                    iy = int(y)
                    if 0 <= ix < width and 0 <= iy < height:
                        mask[iy, ix] = 1

        # Cache mask in metadata and return
        self.ms.meta.mask = mask
    
    @abstractmethod
    def pre_load_meta(self,*args, **kwargs):
        """
        Pre-load metadata before loading full data.

        This method should be called before `load_full_data_from_file` to ensure
        that the necessary metadata is available for data loading.

        Args:
            None

        Returns:
            None
        """
        pass

    @abstractmethod
    def loading_meta(self,*args, **kwargs):
        """
        Pre-load metadata before loading full data.

        This method should be called before `load_full_data_from_file` to ensure
        that the necessary metadata is available for data loading.

        Args:
            None

        Returns:
            None
        """
        pass

    def loaded_meta(self,*args, **kwargs):
        """
        Pre-load metadata before loading full data.

        This method should be called before `load_full_data_from_file` to ensure
        that the necessary metadata is available for data loading.

        Args:
            None

        Returns:
            None
        """
        logger.info("creating ms mask.")
        self.create_ms_meta_mask()
