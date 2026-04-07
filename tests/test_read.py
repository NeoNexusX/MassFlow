import time
import os
import subprocess
from itertools import zip_longest

import pytest
import numpy as np
from massflow.data_manager import MSDataManagerImzML
from massflow.tools.logger import get_logger

logger = get_logger("test_read")

ROUND = 1

READ_BENCHMARK_CASES = [
    # batch_size = 32
    ("batch", "data/example.imzML", 32, 2),
    ("flat", "data/example.imzML", 32, 2),
    ("flat", "data/example.imzML", 32, 8),
    ("matrix", "data/example.imzML", 32, 2),
    ("matrix", "data/example.imzML", 32, 8),
    # batch_size = 64
    ("batch", "data/example.imzML", 64, 2),
    ("batch", "data/example.imzML", 64, 8),
    ("flat", "data/example.imzML", 64, 2),
    ("flat", "data/example.imzML", 64, 8),
    ("matrix", "data/example.imzML", 64, 2),
    ("matrix", "data/example.imzML", 64, 8),
    # batch_size = 128
    ("batch", "data/example.imzML", 128, 2),
    ("batch", "data/example.imzML", 128, 8),
    ("flat", "data/example.imzML", 128, 2),
    ("flat", "data/example.imzML", 128, 8),
    ("matrix", "data/example.imzML", 128, 2),
    ("matrix", "data/example.imzML", 128, 8),
    # batch_size = 256
    ("batch", "data/example.imzML", 256, 2),
    ("batch", "data/example.imzML", 256, 8),
    ("flat", "data/example.imzML", 256, 2),
    ("flat", "data/example.imzML", 256, 8),
    ("matrix", "data/example.imzML", 256, 2),
    ("matrix", "data/example.imzML", 256, 8),
]

def drop_caches():
    if os.name == "nt":
        return
    subprocess.run(["sudo", "purge"], check=True)

def data_reading_benchmark(data_file_path="data/example.imzML"):
    with MSDataManagerImzML(filepath=data_file_path) as ms_md:
        ms_md.load_head_data()
        for spectrum_data in ms_md.ms:
            assert spectrum_data.intensity is not None


def data_reading_full_benchmark(
    read_mode="batch", data_file_path="data/example.imzML", batch_size=1024, max_workers=2
):
    with MSDataManagerImzML(filepath=data_file_path) as ms_md:
        ms_md.load_head_data()

        if read_mode == "batch":
            batch_generator = ms_md.batch_generator(
                batch_size=batch_size,
                max_threads=max_workers
            )
            for batch in batch_generator:
                assert len(batch[1].mz_list) > 0  # type: ignore
            drop_caches()
            return

        if read_mode == "matrix":
            matrix_generator = ms_md.matrix_generator(
                batch_size=batch_size,
                include_mz=True,
                max_threads=max_workers,
            )
            for _mz_data, intensity_matrix, lengths, _coords in matrix_generator:
                assert len(intensity_matrix) > 0
                assert len(lengths) > 0
            drop_caches()
            return

        if read_mode == "flat":
            flat_generator = ms_md.flat_generator(
                batch_size=batch_size,
                include_mz=True,
                max_threads=max_workers
            )
            for _mz_data, intensity_flat, lengths, _coords in flat_generator:
                assert len(intensity_flat) > 0
                assert len(lengths) > 0
            drop_caches()
            return

        raise ValueError(f"Unsupported read_mode: {read_mode}")


def data_reading_by_mode_benchmark(
    read_mode="batch", data_file_path="data/example.imzML", batch_size=1024, max_workers=2
):
    data_reading_full_benchmark(read_mode, data_file_path, batch_size, max_workers)


def _normalize_coordinate(coord) -> tuple:
    if hasattr(coord, "get_tuple"):
        return coord.get_tuple()
    if hasattr(coord, "tolist"):
        return tuple(coord.tolist())
    return tuple(coord)


def _iter_matrix_spectra(matrix_batch, is_continuous: bool):
    mz_data, intensity_matrix, lengths, coords = matrix_batch
    lengths_array = np.asarray(lengths, dtype=np.int64)
    for index, length in enumerate(lengths_array):
        spectrum_length = int(length)
        coord = _normalize_coordinate(coords[index])
        mz_array = np.asarray(mz_data[:spectrum_length]) if is_continuous else np.asarray(mz_data[index, :spectrum_length])
        intensity_array = np.asarray(intensity_matrix[index, :spectrum_length])
        yield coord, spectrum_length, mz_array, intensity_array


def _iter_flat_spectra(flat_batch, is_continuous: bool):
    mz_data, intensity_flat, lengths, coords = flat_batch
    lengths_array = np.asarray(lengths, dtype=np.int64)
    offsets = np.zeros(len(lengths_array), dtype=np.int64)
    if len(lengths_array) > 1:
        offsets[1:] = np.cumsum(lengths_array[:-1], dtype=np.int64)

    for index, length in enumerate(lengths_array):
        start = int(offsets[index])
        end = start + int(length)
        coord = _normalize_coordinate(coords[index])
        mz_array = np.asarray(mz_data[:int(length)]) if is_continuous else np.asarray(mz_data[start:end])
        intensity_array = np.asarray(intensity_flat[start:end])
        yield coord, int(length), mz_array, intensity_array


def _stream_compare_read_modes(data_file_path: str, batch_size: int, max_workers: int) -> None:
    with (
        MSDataManagerImzML(filepath=data_file_path) as batch_dm,
        MSDataManagerImzML(filepath=data_file_path) as matrix_dm,
        MSDataManagerImzML(filepath=data_file_path) as flat_dm,
    ):
        batch_dm.load_head_data()
        matrix_dm.load_head_data()
        flat_dm.load_head_data()

        is_continuous = bool(batch_dm.ms.meta and batch_dm.ms.meta.continuous)

        batch_iter = batch_dm.batch_generator(batch_size=batch_size, max_threads=max_workers)
        matrix_iter = matrix_dm.matrix_generator(
            batch_size=batch_size,
            include_mz=True,
            max_threads=max_workers,
        )
        flat_iter = flat_dm.flat_generator(
            batch_size=batch_size,
            include_mz=True,
            max_threads=max_workers,
        )

        for batch_batch, matrix_batch, flat_batch in zip_longest(batch_iter, matrix_iter, flat_iter):
            assert batch_batch is not None
            assert matrix_batch is not None
            assert flat_batch is not None

            matrix_pixels = list(_iter_matrix_spectra(matrix_batch, is_continuous))
            flat_pixels = list(_iter_flat_spectra(flat_batch, is_continuous))

            assert len(batch_batch) == len(matrix_pixels) == len(flat_pixels)

            for spectrum, matrix_pixel, flat_pixel in zip(batch_batch, matrix_pixels, flat_pixels):
                batch_coord = spectrum.coordinate.get_tuple()
                batch_mz = np.asarray(spectrum.mz_list)
                batch_intensity = np.asarray(spectrum.intensity)

                matrix_coord, matrix_length, matrix_mz, matrix_intensity = matrix_pixel
                flat_coord, flat_length, flat_mz, flat_intensity = flat_pixel

                assert batch_coord == matrix_coord == flat_coord

                valid_length = min(len(batch_mz), len(batch_intensity), matrix_length, flat_length)
                assert len(batch_mz) >= valid_length
                assert len(batch_intensity) >= valid_length

                np.testing.assert_allclose(batch_mz[:valid_length], matrix_mz[:valid_length], rtol=1e-5, atol=1e-8)
                np.testing.assert_allclose(batch_intensity[:valid_length], matrix_intensity[:valid_length], rtol=1e-5, atol=1e-8)
                np.testing.assert_allclose(batch_mz[:valid_length], flat_mz[:valid_length], rtol=1e-5, atol=1e-8)
                np.testing.assert_allclose(batch_intensity[:valid_length], flat_intensity[:valid_length], rtol=1e-5, atol=1e-8)


class TestDataReading:

    """
    Test suite for benchmarking data reading performance of MSDataManagerImzML.
    """

    @pytest.mark.parametrize("data_file_path", ["./data/example.imzML"])
    @pytest.mark.benchmark(timer=time.perf_counter)
    def test_data_reading_benchmark(self, benchmark, data_file_path):
        logger.info("Starting data reading benchmark")
        benchmark.pedantic(
            data_reading_benchmark,
            args=(data_file_path,),
            rounds=ROUND,
            iterations=1,
            warmup_rounds=1,
        )
        logger.info("Finished data reading benchmark")

    @pytest.mark.parametrize(
        "read_mode,data_file_path,batch_size,max_workers",
        READ_BENCHMARK_CASES,
    )
    @pytest.mark.benchmark(timer=time.perf_counter)
    def test_multi_reading_benchmark(self, benchmark, read_mode, data_file_path, batch_size, max_workers):
        logger.info(
            "Starting %s reading benchmark (batch_size=%s, max_workers=%s)",
            read_mode,
            batch_size,
            max_workers,
        )
        benchmark.pedantic(
            data_reading_by_mode_benchmark,
            args=(read_mode, data_file_path, batch_size, max_workers),
            rounds=ROUND,
            iterations=1,
            warmup_rounds=1,
        )
        logger.info(
            "Finished %s reading benchmark (batch_size=%s, max_workers=%s)",
            read_mode,
            batch_size,
            max_workers,
        )

    def test_reading_modes_consistent_by_coordinate(self):
        data_file_path = "data/example.imzML"
        batch_size = 64
        max_workers = 2
        _stream_compare_read_modes(data_file_path, batch_size, max_workers)
