import time
import os
from functools import partial
import pytest
import numpy as np
from massflow.data_manager import MSDataManagerImzML
from massflow.tools.logger import get_logger

logger = get_logger("test_swap")

ROUNDS = 5


def data_swap_out_direct(raw_dm):
    """
    Test: In-memory data -> Write to disk (Swap Out), using writer directly
    """

    swap_dm = MSDataManagerImzML(temp_dir="./temp_swap")
    swap_dm.copy_meta(raw_dm)
    writer = swap_dm.writer

    for raw_spec in raw_dm.ms:

        writer.add_spectrum(
            raw_spec.mz_list, raw_spec.intensity, raw_spec.coordinate.get_tuple()
        )

    swap_dm.close_writer()
    return swap_dm


def cache_matrix_batches(raw_dm, batch_size=512, max_threads=8):
    """Read matrix batches once and keep them in memory for swap-out benchmarks."""
    cached = []
    matrix_generator = raw_dm.matrix_generator(
        batch_size=batch_size,
        include_mz=True,
        max_threads=max_threads,
    )

    for mz_data, intensity_matrix, lengths, coordinates in matrix_generator:
        mz_cache = None if mz_data is None else np.asarray(mz_data).copy()
        cached.append(
            (
                mz_cache,
                np.asarray(intensity_matrix).copy(),
                np.asarray(lengths).copy(),
                np.asarray(coordinates).copy(),
            )
        )

    return cached


def cache_flat_batches(raw_dm, batch_size=512, max_threads=8):
    """Read flat batches once and keep them in memory for swap-out benchmarks."""
    cached = []
    flat_generator = raw_dm.flat_generator(
        batch_size=batch_size,
        include_mz=True,
        max_threads=max_threads,
    )

    for mz_data, intensity_flat, lengths, coordinates in flat_generator:
        mz_cache = None if mz_data is None else np.asarray(mz_data).copy()
        cached.append(
            (
                mz_cache,
                np.asarray(intensity_flat).copy(),
                np.asarray(lengths).copy(),
                np.asarray(coordinates).copy(),
            )
        )

    return cached


def data_swap_out_matrix_cached(raw_dm, cached_matrix_batches):
    """Swap out using pre-cached matrix batches to isolate write performance."""
    swap_dm = MSDataManagerImzML(temp_dir="./temp_swap")
    swap_dm.copy_meta(raw_dm)

    for mz_data, intensity_matrix, lengths, coordinates in cached_matrix_batches:
        mz_batch = mz_data

        swap_dm.swap_matrix_data_out2disk(
            mz_data=mz_batch,
            intensity_matrix=intensity_matrix,
            lengths=lengths,
            coordinates=coordinates,
        )

    swap_dm.close_writer()
    return swap_dm


def data_swap_out_flat_cached(raw_dm, cached_flat_batches):
    """Swap out using pre-cached flat batches to isolate write performance."""
    is_continuous = bool(raw_dm.ms.meta and raw_dm.ms.meta.continuous)

    swap_dm = MSDataManagerImzML(temp_dir="./temp_swap")
    swap_dm.copy_meta(raw_dm)

    for mz_data, intensity_flat, lengths, coordinates in cached_flat_batches:
        lengths_batch = lengths
        point_count = int(np.sum(lengths_batch, dtype=np.int64))
        intensity_batch = intensity_flat[:point_count]

        if mz_data is None:
            mz_batch = None
        elif is_continuous:
            mz_batch = mz_data
        else:
            mz_batch = mz_data[:point_count]

        swap_dm.swap_flat_data_out2disk(
            mz_flat=mz_batch,
            intensity_flat=intensity_batch,
            lengths=lengths_batch,
            coordinates=coordinates,
        )

    swap_dm.close_writer()
    return swap_dm


def iter_flat_cached_spectra(cached_flat_batches, is_continuous):
    """Yield per-spectrum tuples (coord, mz, intensity) from cached flat batches."""
    for mz_data, intensity_flat, lengths, coordinates in cached_flat_batches:
        lengths_array = np.asarray(lengths, dtype=np.int64)
        offsets = np.zeros(len(lengths_array), dtype=np.int64)
        if len(lengths_array) > 1:
            offsets[1:] = np.cumsum(lengths_array[:-1], dtype=np.int64)

        for index, length in enumerate(lengths_array):
            start = int(offsets[index])
            end = start + int(length)
            coord = tuple(int(v) for v in np.asarray(coordinates[index]).tolist())

            if mz_data is None:
                mz_array = None
            elif is_continuous:
                mz_array = np.asarray(mz_data[: int(length)])
            else:
                mz_array = np.asarray(mz_data[start:end])

            intensity_array = np.asarray(intensity_flat[start:end])
            yield coord, mz_array, intensity_array


def data_swap_in_benchmark(swap_dm, mode="multithread"):
    """
    Test: Disk file -> Load and access in memory (Swap In)
    """
    dm = swap_dm
    dm.load_full_data_from_file()
    dm.inspect_data()

    if mode == "multithread":
        batch_generator = dm.batch_generator(batch_size=512, max_threads=32)

        for batch in batch_generator:
            assert len(batch) > 0

            # clean batch memory
            for spec in batch:
                spec.clear_data()
    else:

        # read in memory
        for spec in dm.ms:
            assert spec.intensity is not None

        # clean memory
        for spec in dm.ms:
            spec.clear_data()

    logger.info("Completed loading and accessing swapped data.")


class TestDataSwap:
    """
    Disk swap performance test class
    """

    @pytest.fixture(scope="class")
    def ms_raw_data(self):

        filepath = r"./data/example.imzML"

        if not os.path.exists(filepath):
            pytest.skip(f"Original data not found at {filepath}")

        # load meta data
        dm_raw = MSDataManagerImzML(filepath=filepath)
        dm_raw.load_head_data()
        dm_raw.inspect_data()
        # load data check
        batch_generator = dm_raw.batch_generator()
        for batch in batch_generator:
            assert len(batch) > 0

        return dm_raw

    @pytest.fixture(scope="class")
    def matrix_cached_batches(self, ms_raw_data):
        cached = cache_matrix_batches(ms_raw_data)
        assert len(cached) > 0
        return cached

    @pytest.fixture(scope="class")
    def flat_cached_batches(self, ms_raw_data):
        cached = cache_flat_batches(ms_raw_data)
        assert len(cached) > 0
        return cached

    def test_swap_out_in_data_consistency_flat_cached(
        self, ms_raw_data, flat_cached_batches
    ):
        """Validate flat cached swap-out then swap-in preserves coord/mz/intensity."""
        is_continuous = bool(ms_raw_data.ms.meta and ms_raw_data.ms.meta.continuous)

        expected_by_coord = {}
        for coord, mz_array, intensity_array in iter_flat_cached_spectra(
            flat_cached_batches, is_continuous
        ):
            expected_by_coord[coord] = (
                None if mz_array is None else np.asarray(mz_array).copy(),
                np.asarray(intensity_array).copy(),
            )

        assert len(expected_by_coord) > 0

        dm_swapped = data_swap_out_flat_cached(ms_raw_data, flat_cached_batches)
        try:
            dm_swapped.load_head_data()
            logger.info(f"{dm_swapped.file_base_path} loaded for swap-in consistency check")
            logger.info(f"Original spectra count: {len(expected_by_coord)}, Swapped spectra count: {len(dm_swapped.ms)}")

            for spectrum in dm_swapped.ms:
                coord = spectrum.coordinate.get_tuple()
                assert coord in expected_by_coord
                expected_mz, expected_intensity = expected_by_coord.pop(coord)
                actual_mz = None if spectrum.mz_list is None else np.asarray(spectrum.mz_list)
                actual_intensity = np.asarray(spectrum.intensity)

                np.testing.assert_allclose(
                    actual_intensity,
                    expected_intensity,
                    rtol=1e-5,
                    atol=1e-8,
                )

                if expected_mz is None:
                    assert actual_mz is None
                else:
                    assert actual_mz is not None
                    np.testing.assert_allclose(
                        actual_mz,
                        expected_mz,
                        rtol=1e-5,
                        atol=1e-8,
                    )

            assert len(expected_by_coord) == 0
        finally:
            dm_swapped.close()

    @pytest.mark.benchmark(timer=time.perf_counter)
    def test_swap_out_performance_direct(self, benchmark, ms_raw_data):
        """
        Performance test: Data write speed to disk
        """
        raw_dm = ms_raw_data

        logger.info("Starting Swap Out benchmark for full dataset")

        dm_result = benchmark.pedantic(
            data_swap_out_direct,
            args=(raw_dm,),
            rounds=ROUNDS,
            iterations=1,
            warmup_rounds=1,
        )

        dm_result.close()
        logger.info("Finished Swap Out benchmark")

    @pytest.mark.benchmark(timer=time.perf_counter)
    def test_swap_out_performance_matrix_cached(
        self, benchmark, ms_raw_data, matrix_cached_batches
    ):
        """Performance test: matrix read cache -> swap out write speed."""
        logger.info("Starting Matrix Cached Swap Out benchmark for full dataset")

        dm_result = benchmark.pedantic(
            data_swap_out_matrix_cached,
            args=(ms_raw_data, matrix_cached_batches),
            rounds=ROUNDS,
            iterations=1,
            warmup_rounds=1,
        )

        dm_result.close()
        logger.info("Finished Matrix Cached Swap Out benchmark")

    @pytest.mark.benchmark(timer=time.perf_counter)
    def test_swap_out_performance_flat_cached(
        self, benchmark, ms_raw_data, flat_cached_batches
    ):
        """Performance test: flat read cache -> swap out write speed."""
        logger.info("Starting Flat Cached Swap Out benchmark for full dataset")

        dm_result = benchmark.pedantic(
            data_swap_out_flat_cached,
            args=(ms_raw_data, flat_cached_batches),
            rounds=ROUNDS,
            iterations=1,
            warmup_rounds=1,
        )

        dm_result.close()
        logger.info("Finished Flat Cached Swap Out benchmark")


    @pytest.mark.benchmark(timer=time.perf_counter)
    @pytest.mark.parametrize("mode", ["normal", "multithread"])
    def test_swap_in_performance(self, benchmark, ms_raw_data, mode):
        """
        Performance test: Data read speed from disk
        """
        logger.info("Starting Swap In benchmark for full dataset")

        raw_dm = ms_raw_data

        def set_up_swap_data():
            dm_swap_new = data_swap_out_direct(raw_dm)
            raw_dm.clear_all_data_memory()
            return (dm_swap_new,), {}

        benchmark.pedantic(
            partial(data_swap_in_benchmark, mode=mode),
            setup=set_up_swap_data,
            rounds=ROUNDS,
            iterations=1,
            warmup_rounds=1,
        )

        logger.info("Finished Swap In benchmark")
