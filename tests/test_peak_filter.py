import numpy as np

from massflow.preprocess.helper.peak_align_helper import (
    _filter_reference_by_counts,
    bin_peaks_parallel,
    merge_peaks,
)
from massflow.preprocess.helper.peak_filter_helper import (
    PeakReferenceStats,
    filter_reference_peaks,
)


class _FakePeakManager:
    def __init__(self):
        self.ms = [object(), object(), object()]

    def flat_generator(self, **_kwargs):
        mz_flat = np.array([100.0, 101.0, 100.0, 102.0, 100.0], dtype=np.float64)
        intensity = np.ones(mz_flat.size, dtype=np.float64)
        lengths = np.array([2, 2, 1], dtype=np.int32)
        coordinates = np.zeros((3, 3), dtype=np.int32)
        yield mz_flat, intensity, lengths, coordinates


def test_bin_peaks_parallel_returns_reference_counts():
    reference, counts = bin_peaks_parallel(
        data_manager=_FakePeakManager(), # type: ignore
        domain=np.array([100.0, 101.0, 102.0], dtype=np.float64),
        tolerance=0.01,
        tol_method="abs",
        batch_size=2,
        max_threads=0,
    )

    np.testing.assert_array_equal(reference, np.array([100.0, 101.0, 102.0]))
    np.testing.assert_array_equal(counts, np.array([3, 1, 1]))


def test_merge_peaks_returns_counts_aligned_with_reference():
    mz_values = np.array([100.000, 100.005, 100.010, 101.000], dtype=np.float64)
    counts = np.array([2, 5, 1, 3], dtype=np.int64)

    merged_mz, merged_counts = merge_peaks(
        mz_list=mz_values.copy(),
        counts=counts.copy(),
        tolerance=0.02,
        tol_method="abs",
    )

    expected_first = (100.000 * 2 + 100.005 * 5 + 100.010 * 1) / 8
    np.testing.assert_allclose(merged_mz, np.array([expected_first, 101.000]))
    np.testing.assert_array_equal(merged_counts, np.array([8, 3]))


def test_filter_reference_by_counts_supports_count_and_frequency():
    reference = np.array([100.0, 101.0, 102.0, 103.0])
    counts = np.array([1, 5, 10, 20], dtype=np.int64)

    filtered, filtered_counts = _filter_reference_by_counts(
        reference,
        counts,
        n_spectra=20,
        min_count=5,
        min_frequency=0.5,
    )

    np.testing.assert_array_equal(filtered, np.array([102.0, 103.0]))
    np.testing.assert_array_equal(filtered_counts, np.array([10, 20]))


def test_frequency_is_capped_when_merged_counts_exceed_spectrum_count():
    stats = PeakReferenceStats(
        reference=np.array([100.0]),
        counts=np.array([15], dtype=np.int64),
        n_spectra=10,
    )

    np.testing.assert_array_equal(stats.frequency, np.array([1.0]))


def test_filter_reference_peaks_combines_frequency_and_intensity_rules():
    stats = PeakReferenceStats(
        reference=np.array([100.0, 101.0, 102.0, 103.0]),
        counts=np.array([2, 8, 10, 1], dtype=np.int64),
        intensity_sum=np.array([10.0, 20.0, 80.0, 100.0]),
        intensity_max=np.array([6.0, 9.0, 30.0, 100.0]),
        n_spectra=10,
    )

    filtered, keep = filter_reference_peaks(
        stats,
        min_frequency=0.5,
        min_mean_intensity=5.0,
        min_max_intensity=20.0,
    )

    np.testing.assert_array_equal(filtered, np.array([102.0]))
    np.testing.assert_array_equal(keep, np.array([False, False, True, False]))


def test_filter_reference_peaks_top_k_preserves_original_order():
    stats = PeakReferenceStats(
        reference=np.array([100.0, 101.0, 102.0, 103.0]),
        counts=np.array([1, 5, 3, 2], dtype=np.int64),
    )

    filtered, keep = filter_reference_peaks(stats, top_k=2, top_k_by="count")

    np.testing.assert_array_equal(filtered, np.array([101.0, 102.0]))
    np.testing.assert_array_equal(keep, np.array([False, True, True, False]))


def test_filter_reference_peaks_can_remove_simple_isotope_candidates():
    stats = PeakReferenceStats(
        reference=np.array([500.0, 501.0033548378, 502.5]),
        counts=np.array([10, 4, 8], dtype=np.int64),
    )

    filtered, keep = filter_reference_peaks(
        stats,
        remove_isotopes=True,
        isotope_score="count",
        isotope_charge_states=(1,),
        isotope_mz_tolerance=0.001,
    )

    np.testing.assert_array_equal(filtered, np.array([500.0, 502.5]))
    np.testing.assert_array_equal(keep, np.array([True, False, True]))
