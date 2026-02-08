import time
import os
import pytest
import numpy as np
from massflow.module.mass_spectrum_set import MassSpectrumSet
from massflow.module.ms_data_manager_imzml import MSDataManagerImzML
from massflow.tools.logger import get_logger

logger = get_logger("test_matrix_read_swap")

ROUND = 1

def data_read2matrix(dm, batch_size, include_mz, max_threads, swap):
    """
    Test: Disk data -> Read to matrix, with/without swap.
    """
    swap_dm = MSDataManagerImzML(MassSpectrumSet(), temp_dir="./tmp")
    swap_dm.copy_meta(dm)

    for mz_data, intensity_matrix, batch_lengths, batch_coordinates in dm.get_matrix_generator(
        batch_size=batch_size,
        include_mz=include_mz,
        max_threads=max_threads,
    ):

        if swap:
            swap_dm.swap_matrix_data_out2disk(mz_data, intensity_matrix, batch_lengths, batch_coordinates)
        else:
            _ = intensity_matrix.sum()

    if swap:
        swap_dm.close_writer()
        return swap_dm

def data_read2batch(dm, batch_size, max_threads, swap):
    """
    Test: Disk data -> Read to batch, with/without swap.
    """
    swap_dm = MSDataManagerImzML(MassSpectrumSet(), temp_dir="./tmp")
    swap_dm.copy_meta(dm)
    writer = swap_dm.writer

    for batch in dm.get_batch_generator(
        batch_size=batch_size,
        max_threads=max_threads,
    ):

        if swap:
            for spec in batch:
                writer.addSpectrum(spec.mz_list, spec.intensity, spec.coordinate.get_tuple())
        else:
            _ = batch

        dm.clear_batch_data_memory(batch)

    if swap:
        swap_dm.close_writer()
        return swap_dm

def _collect_batch_data(dm):
    """
    collect all data via get_batch_generator.
    """
    results = []
    for batch in dm.get_batch_generator(batch_size=1024, max_threads=8):
        for spec in batch:
            results.append({
                "mz": np.array(spec.mz_list, copy=True),
                "intensity": np.array(spec.intensity, copy=True),
                "coord": spec.coordinate.get_tuple(),
            })
        dm.clear_batch_data_memory(batch)
    return results

def _collect_matrix_data(dm):
    """
    collect all data via get_matrix_generator (include_mz=True).
    """
    results = []
    for mz_data, intensity_matrix, lengths, coordinates in dm.get_matrix_generator(
        batch_size=1024, include_mz=True, max_threads=8,
    ):
        n = intensity_matrix.shape[0]
        is_continuous = (mz_data is not None and mz_data.ndim == 1)
        for i in range(n):
            length = int(lengths[i])
            intensity = intensity_matrix[i, :length].copy()
            mz = mz_data[:length].copy() if is_continuous else mz_data[i, :length].copy()
            coord = tuple(int(c) for c in coordinates[i])

            results.append({
                "mz": mz,
                "intensity": intensity,
                "coord": coord,
            })
    return results

def _compare_spectrum_lists(list_a, list_b):
    """
    Strictly compare the consistency of two lists of spectrum data.

    Check items:
        1. Total number of spectra matches.
        2. Coordinates of each spectrum match.
        3. m/z array length of each spectrum matches.
        4. Intensity array length of each spectrum matches.
        5. m/z values match element-wise (floating point tolerance atol=1e-6, rtol=1e-6).
        6. Intensity values match element-wise (floating point tolerance atol=1e-4, rtol=1e-5).

    Args:
        list_a: First set of data.
        list_b: Second set of data.

    Raises:
        AssertionError: If any check fails.
    """
    # Total number matches
    assert len(list_a) == len(list_b), (f"Spectrum count mismatch: batch_read={len(list_a)}, matrix_read={len(list_b)}")

    # Compare spectra one by one in order
    for idx, (a, b) in enumerate(zip(list_a, list_b)):

        # Check coordinate consistency; m/z and intensity lengths match
        assert a["coord"] == b["coord"], (f"spectrum #{idx} (batch_read vs matrix_read): Coordinate mismatch {a['coord']} != {b['coord']}")
        assert len(a["mz"]) == len(b["mz"]), (f"spectrum #{idx} (batch_read vs matrix_read): m/z length mismatch {len(a['mz'])} != {len(b['mz'])}")
        assert len(a["intensity"]) == len(b["intensity"]), (
            f"spectrum #{idx} (batch_read vs matrix_read): intensity length mismatch {len(a['intensity'])} != {len(b['intensity'])}"
        )

        # Compare m/z values
        np.testing.assert_allclose(
            a["mz"], b["mz"],
            atol=1e-6, rtol=1e-6,
            err_msg=f"spectrum #{idx} (batch_read vs matrix_read): m/z value mismatch"
        )

        # Compare intensity values
        np.testing.assert_allclose(
            a["intensity"], b["intensity"],
            atol=1e-4, rtol=1e-5,
            err_msg=f"spectrum #{idx} (batch_read vs matrix_read): intensity value mismatch"
        )

    logger.info(f"Data consistency verification passed (batch_read vs matrix_read): All {len(list_a)} spectra matched.")

def _compare_by_coordinate_map(list_a, list_b):
    """
    Compare two sets of spectrum data based on coordinate mapping (order-independent).

    Check items:
        1. Total number of spectra matches.
        2. Coordinate sets match completely (no missing, no extras).
        3. m/z values at the same coordinate match.
        4. Intensity values at the same coordinate match.

    Raises:
        AssertionError: If any check fails.
    """
    assert len(list_a) == len(list_b), (
        f"Spectrum count mismatch: original={len(list_a)}, matrix_swap_readback={len(list_b)}"
    )

    # Build coordinate -> data mapping
    map_a = {item["coord"]: item for item in list_a}
    map_b = {item["coord"]: item for item in list_b}

    # Check consistency of coordinate sets
    coords_a = set(map_a.keys())
    coords_b = set(map_b.keys())
    only_in_a = coords_a - coords_b
    only_in_b = coords_b - coords_a
    assert coords_a == coords_b, (
        f"Coordinate set mismatch: "
        f"Only in original: {list(only_in_a)[:5]}{'...' if len(only_in_a) > 5 else ''}, "
        f"Only in matrix_swap_readback: {list(only_in_b)[:5]}{'...' if len(only_in_b) > 5 else ''}"
    )

    # Compare coordinate by coordinate
    for coord in sorted(coords_a):
        a = map_a[coord]
        b = map_b[coord]

        # Length
        assert len(a["mz"]) == len(b["mz"]), (
            f"coord {coord} (original vs matrix_swap_readback): m/z length mismatch {len(a['mz'])} != {len(b['mz'])}"
        )
        assert len(a["intensity"]) == len(b["intensity"]), (
            f"coord {coord} (original vs matrix_swap_readback): intensity length mismatch {len(a['intensity'])} != {len(b['intensity'])}"
        )

        # m/z and intensity values
        np.testing.assert_allclose(
            a["mz"], b["mz"],
            atol=1e-6, rtol=1e-6,
            err_msg=f"coord {coord} (original vs matrix_swap_readback): m/z value mismatch"
        )
        np.testing.assert_allclose(
            a["intensity"], b["intensity"],
            atol=1e-4, rtol=1e-5,
            err_msg=f"coord {coord} (original vs matrix_swap_readback): intensity value mismatch"
        )

    logger.info(
        f"Data consistency verification based on coordinate mapping passed (original vs matrix_swap_readback): "
        f"All {len(list_a)} spectra, {len(coords_a)} unique coordinates matched."
    )

def verify_read_data_consistency(dm):
    """
    Verify data consistency between batch read mode and matrix read mode.

    Core logic:
        1. Collect all mz, intensity, coordinate spectrum by spectrum via get_batch_generator.
        2. Collect the same data row by row via get_matrix_generator (include_mz=True).
        3. Strictly compare the two for complete consistency in spectrum count, order, coordinates, m/z, and intensity.

    This is the basic verification of data correctness: ensuring that both read channels produce exactly the same results for the same original file.
    """
    # Collect batch data
    batch_data = _collect_batch_data(dm)
    matrix_data = _collect_matrix_data(dm)

    _compare_spectrum_lists(batch_data, matrix_data)
    logger.info("READ consistency verification completed: batch and matrix read data are completely consistent.")

def verify_swap_data_consistency(dm):
    """
    Verify data correctness after swap (write to disk and read back).

    Core logic:
        1. Collect baseline data (ground truth) from original data via the batch path.
        2. Swap out to disk via the matrix path -> read back -> collect data via the batch path.
        3. Compare: Original vs matrix-swap read back -> Verify data fidelity of matrix swap write/read back.

    Use coordinate mapping for comparison because the order of spectra may differ when reading back after swap.
    """
    logger.info("Start verifying SWAP data consistency (original vs swap read back)")

    # Collect original baseline data
    original_data = _collect_batch_data(dm)

    # matrix swap -> read back
    matrix_swap_dm = data_read2matrix(dm, batch_size=1024, include_mz=True, max_threads=8, swap=True)
    matrix_swap_dm.load_full_data_from_file()
    matrix_swap_readback = _collect_batch_data(matrix_swap_dm)

    # Original vs matrix-swap read back
    _compare_by_coordinate_map(original_data, matrix_swap_readback)

    # Clean up temporary resources
    matrix_swap_dm.close()

    logger.info("SWAP consistency verification completed: matrix-swap read back is completely consistent with original data, test passed.")

class TestMatrixGenerator:
    """
    Test the get_matrix_generator method of MSDataManagerImzML.
    """

    @pytest.fixture(scope="class")
    def ms_raw_data(self):

        filepath = r"E:\data\data\example\example.imzML"
        if not os.path.exists(filepath):
            pytest.skip(f"Original data not found at {filepath}")
        # loda meta data
        dm_raw = MSDataManagerImzML(filepath=filepath)
        dm_raw.load_full_data_from_file()

        return dm_raw

    @pytest.mark.parametrize("batch_size", [1024])
    @pytest.mark.parametrize("max_threads", [8])
    @pytest.mark.parametrize(("include_mz", "swap"), [(True, False), (True, True)]) # swap=True need include_mz=True
    @pytest.mark.benchmark(timer=time.perf_counter)
    def test_matrix_read_swap_performance(self, benchmark, ms_raw_data, batch_size, max_threads, include_mz, swap):
        """
        Performance test: Data read to matrix with/without swap.
        """

        dm_result = benchmark.pedantic(
            data_read2matrix,
            args=(ms_raw_data, batch_size, include_mz, max_threads, swap),
            rounds=ROUND,
            iterations=1,
            warmup_rounds=1
        )

        if swap:
            dm_result.close()

    @pytest.mark.parametrize("batch_size", [1024])
    @pytest.mark.parametrize("max_threads", [8])
    @pytest.mark.parametrize("swap", [True, False])
    @pytest.mark.benchmark(timer=time.perf_counter)
    def test_batch_read_swap_performance(self, benchmark, ms_raw_data, batch_size, max_threads, swap):
        """
        Performance test: Data read to batch with/without swap.
        """

        dm_result = benchmark.pedantic(
            data_read2batch,
            args=(ms_raw_data, batch_size, max_threads, swap),
            rounds=ROUND,
            iterations=1,
            warmup_rounds=1
        )

        if swap:
            dm_result.close()

    # 如果数据比较大的话，下面的测试会比较耗时，且占用内存比较大

    # def test_verify_read_data_consistency(self, ms_raw_data):
    #     """
    #     Verify data consistency between batch and matrix read (without swap).
    #     """

    #     verify_read_data_consistency(ms_raw_data)

    # def test_verify_swap_data_consistency(self, ms_raw_data):
    #     """
    #     Verify data correctness after swap (write to disk and read back).
    #     """

    #     verify_swap_data_consistency(ms_raw_data)
