import time
import os
import pytest
import numpy as np
from massflow.module import MassSpectrumSet
from massflow.data_manager import MSDataManagerImzML
from massflow.tools import get_logger

logger = get_logger("test_matrix_read_swap")

ROUND = 2

def data_read2matrix(dm, batch_size, max_threads, swap):
    """Test: Disk data -> Read to matrix, with/without swap."""
    if swap:
        swap_dm = MSDataManagerImzML(MassSpectrumSet(), temp_dir="./tmp")
        swap_dm.copy_meta(dm)

    for mz_data, intensity_matrix, batch_lengths, batch_coordinates in dm.matrix_generator(
        batch_size=batch_size,
        include_mz=True,
        max_threads=max_threads,
    ):

        if swap:
            swap_dm.swap_matrix_data_out2disk(mz_data, intensity_matrix, batch_lengths, batch_coordinates)
        else:
            _ = intensity_matrix

    if swap:
        swap_dm.close_writer()
        return swap_dm

def data_read2batch(dm, batch_size, max_threads, swap):
    """Test: Disk data -> Read to batch, with/without swap."""
    if swap:
        swap_dm = MSDataManagerImzML(MassSpectrumSet(), temp_dir="./tmp")
        swap_dm.copy_meta(dm)
        writer = swap_dm.writer

    for batch in dm.batch_generator(
        batch_size=batch_size,
        max_threads=max_threads,
    ):

        if swap:
            for spec in batch:
                writer.add_spectrum(spec.mz_list, spec.intensity, spec.coordinate.get_tuple())
        else:
            _ = batch

        dm.clear_batch_data_memory(batch)

    if swap:
        swap_dm.close_writer()
        return swap_dm

def _matrix_data_generator(dm):
    """Convert matrix blocks into a stream of individual spectrum tuples."""
    for mz_data, intensity_matrix, lengths, coords in dm.matrix_generator(
        batch_size=256, include_mz=True
    ):
        n_spectra = intensity_matrix.shape[0]
        is_shared_mz = (mz_data is not None and mz_data.ndim == 1)

        for i in range(n_spectra):
            length = int(lengths[i])
            intensity = intensity_matrix[i, :length]
            mz = mz_data[:length] if is_shared_mz else mz_data[i, :length]

            coord = tuple(coords[i].tolist())

            yield (coord, mz, intensity)

def _batch_data_generator(dm):
    """Convert batch objects into a stream of individual spectrum tuples."""
    for batch in dm.batch_generator(batch_size=256):
        for spec in batch:
            yield (spec.coordinate.get_tuple(), spec.mz_list, spec.intensity)

        dm.clear_batch_data_memory(batch)

def verify_data_streaming(dm_original, dm_swap=None):
    """Verify data consistency between original and swap (if provided) by streaming through both generators in parallel."""
    gen_a = _batch_data_generator(dm_original)
    gen_b = _batch_data_generator(dm_swap) if dm_swap else _matrix_data_generator(dm_original)

    count = 0
    for count, (data_a, data_b) in enumerate(zip(gen_a, gen_b), 1):
        coord_a, mz_a, int_a = data_a
        coord_b, mz_b, int_b = data_b

        assert coord_a[:2] == coord_b[:2], f"Coordinate mismatch at index {count-1}: {coord_a} vs {coord_b}"
        np.testing.assert_allclose(mz_a, mz_b, rtol=1e-5, atol=1e-8, err_msg=f"m/z mismatch at index {count-1}")
        np.testing.assert_allclose(int_a, int_b, rtol=1e-5, atol=1e-8, err_msg=f"Intensity mismatch at index {count-1}")

    expected = len(dm_swap.ms) if dm_swap else len(dm_original.ms)
    assert count == expected, f"Spectrum count mismatch: processed {count}, expected {expected}"

    logger.info(f"successful, verified {count} spectra.")

class TestMatrixGenerator:
    """Test the matrix_generator method of MSDataManagerImzML."""

    @pytest.fixture(scope="class")
    def ms_raw_data(self):

        filepath = r"./data/example.imzML"
        if not os.path.exists(filepath):
            pytest.skip(f"Original data not found at {filepath}")

        dm_raw = MSDataManagerImzML(filepath=filepath)
        dm_raw.load_head_data()
        dm_raw.inspect_data()

        return dm_raw

    @pytest.mark.parametrize("batch_size", [256])
    @pytest.mark.parametrize("max_threads", [2])
    @pytest.mark.parametrize(("swap"), [True, False])
    @pytest.mark.benchmark(timer=time.perf_counter)
    def test_matrix_read_swap_performance(self, benchmark, ms_raw_data, batch_size, max_threads, swap):
        """Performance test: Data read to matrix with/without swap."""

        dm_result = benchmark.pedantic(
            data_read2matrix,
            args=(ms_raw_data, batch_size, max_threads, swap),
            rounds=ROUND,
            iterations=1,
            warmup_rounds=0
        )

        if swap:
            dm_result.close()

    @pytest.mark.parametrize("batch_size", [256])
    @pytest.mark.parametrize("max_threads", [2])
    @pytest.mark.parametrize("swap", [True, False])
    @pytest.mark.benchmark(timer=time.perf_counter)
    def test_batch_read_swap_performance(self, benchmark, ms_raw_data, batch_size, max_threads, swap):
        """Performance test: Data read to batch with/without swap."""

        dm_result = benchmark.pedantic(
            data_read2batch,
            args=(ms_raw_data, batch_size, max_threads, swap),
            rounds=ROUND,
            iterations=1,
            warmup_rounds=0
        )

        if swap:
            dm_result.close()

    def test_verify_read_data_consistency(self, ms_raw_data):
        """Verify data consistency between batch and matrix read (without swap)."""

        verify_data_streaming(ms_raw_data, None)

    def test_verify_swap_data_consistency(self, ms_raw_data):
        """Verify data correctness after swap (write to disk and read back)."""
        dm_swap = data_read2matrix(ms_raw_data, batch_size=256, max_threads=2, swap=True)
        try:
            dm_swap.load_head_data()
            verify_data_streaming(ms_raw_data, dm_swap)
        finally:
            dm_swap.close()
