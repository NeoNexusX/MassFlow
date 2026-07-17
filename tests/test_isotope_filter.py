"""Algorithm-level tests for isotope peak filtering."""

import numpy as np
import pytest

from massflow.preprocess.helper.isotope_filter_helper import (
    filter_peaks_by_isotopes,
)
from massflow.preprocess.numba.isotope_filter_numba import ISOTOPE_MASS_DIFFERENCE


DELTA = ISOTOPE_MASS_DIFFERENCE


# ---------------------------------------------------------------------------
# Isotope-cluster behavior
# ---------------------------------------------------------------------------


def test_filters_singly_charged_isotope_cluster() -> None:
    mz_data = np.array([100.0, 100.0 + DELTA, 100.0 + 2 * DELTA, 150.0])
    intensity = np.array([100.0, 80.0, 60.0, 5.0], dtype=np.float32)

    filtered_mz, filtered_intensity, filtered_lengths = filter_peaks_by_isotopes(
        mz_data,
        intensity,
        np.array([4]),
        mz_tolerance=1e-5,
        units="da",
    )

    np.testing.assert_allclose(filtered_mz, [100.0, 150.0])
    np.testing.assert_allclose(filtered_intensity, [100.0, 5.0])
    np.testing.assert_array_equal(
        filtered_lengths,
        np.array([2], dtype=np.uint32),
    )
    assert filtered_mz.dtype == mz_data.dtype
    assert filtered_intensity.dtype == intensity.dtype


def test_compacts_each_spectrum_independently() -> None:
    mz_data = np.array(
        [
            100.0,
            100.0 + DELTA,
            100.0 + 2 * DELTA,
            150.0,
            200.0,
            201.2,
            202.4,
        ]
    )
    intensity = np.array([100.0, 80.0, 60.0, 5.0, 30.0, 20.0, 10.0])

    filtered_mz, filtered_intensity, filtered_lengths = filter_peaks_by_isotopes(
        mz_data,
        intensity,
        np.array([4, 3]),
        mz_tolerance=1e-5,
        units="da",
    )

    np.testing.assert_allclose(filtered_mz, [100.0, 150.0, 200.0, 201.2, 202.4])
    np.testing.assert_allclose(filtered_intensity, [100.0, 5.0, 30.0, 20.0, 10.0])
    np.testing.assert_array_equal(filtered_lengths, [2, 3])


@pytest.mark.parametrize(
    ("min_charge", "max_charge", "expected_mz"),
    [
        (1, 1, [500.0, 500.0 + DELTA / 2, 500.0 + DELTA]),
        (2, 2, [500.0]),
    ],
)
def test_respects_charge_range(
    min_charge: int,
    max_charge: int,
    expected_mz: list[float],
) -> None:
    mz_data = np.array([500.0, 500.0 + DELTA / 2, 500.0 + DELTA])
    intensity = np.array([100.0, 80.0, 60.0])

    filtered_mz, _, _ = filter_peaks_by_isotopes(
        mz_data,
        intensity,
        np.array([3]),
        mz_tolerance=1e-5,
        units="da",
        min_charge=min_charge,
        max_charge=max_charge,
    )

    np.testing.assert_allclose(filtered_mz, expected_mz)


def test_selects_nearest_peak_within_tolerance() -> None:
    farther_candidate = 100.0 + DELTA - 0.004
    nearest_candidate = 100.0 + DELTA + 0.001
    mz_data = np.array([100.0, farther_candidate, nearest_candidate, 100.0 + 2 * DELTA])
    intensity = np.array([100.0, 70.0, 80.0, 60.0])

    filtered_mz, filtered_intensity, filtered_lengths = filter_peaks_by_isotopes(
        mz_data,
        intensity,
        np.array([4]),
        mz_tolerance=0.01,
        units="da",
    )

    np.testing.assert_allclose(filtered_mz, [100.0, farther_candidate])
    np.testing.assert_allclose(filtered_intensity, [100.0, 70.0])
    np.testing.assert_array_equal(filtered_lengths, [2])


@pytest.mark.parametrize(
    ("mz_tolerance", "expected_mz"),
    [
        (0.01, [100.0, 100.0 + DELTA + 0.02, 100.0 + 2 * DELTA + 0.02]),
        (0.03, [100.0]),
    ],
)
def test_respects_absolute_tolerance(
    mz_tolerance: float,
    expected_mz: list[float],
) -> None:
    mz_data = np.array([100.0, 100.0 + DELTA + 0.02, 100.0 + 2 * DELTA + 0.02])
    intensity = np.array([100.0, 80.0, 60.0])

    filtered_mz, _, _ = filter_peaks_by_isotopes(
        mz_data,
        intensity,
        np.array([3]),
        mz_tolerance=mz_tolerance,
        units="da",
    )

    np.testing.assert_allclose(filtered_mz, expected_mz)


@pytest.mark.parametrize(
    ("monoisotopic_mz", "expected_length"),
    [(100.0, 3), (1000.0, 1)],
)
def test_ppm_tolerance_scales_with_mz(
    monoisotopic_mz: float,
    expected_length: int,
) -> None:
    mz_data = np.array(
        [
            monoisotopic_mz,
            monoisotopic_mz + DELTA + 0.008,
            monoisotopic_mz + 2 * DELTA + 0.008,
        ]
    )
    intensity = np.array([100.0, 80.0, 60.0])

    filtered_mz, _, filtered_lengths = filter_peaks_by_isotopes(
        mz_data,
        intensity,
        np.array([3]),
        mz_tolerance=10.0,
        units="ppm",
    )

    assert filtered_mz.size == expected_length
    np.testing.assert_array_equal(filtered_lengths, [expected_length])


@pytest.mark.parametrize(
    ("start_intensity_check", "expected_length"),
    [(1, 3), (2, 1)],
)
def test_start_intensity_check_controls_first_compared_isotope(
    start_intensity_check: int,
    expected_length: int,
) -> None:
    mz_data = np.array([100.0, 100.0 + DELTA, 100.0 + 2 * DELTA])
    intensity = np.array([10.0, 20.0, 15.0])

    filtered_mz, _, _ = filter_peaks_by_isotopes(
        mz_data,
        intensity,
        np.array([3]),
        mz_tolerance=1e-5,
        units="da",
        start_intensity_check=start_intensity_check,
    )

    assert filtered_mz.size == expected_length


@pytest.mark.parametrize(
    ("use_decreasing_model", "expected_length"),
    [(True, 3), (False, 1)],
)
def test_decreasing_intensity_model_can_be_disabled(
    use_decreasing_model: bool,
    expected_length: int,
) -> None:
    mz_data = np.array([100.0, 100.0 + DELTA, 100.0 + 2 * DELTA])
    intensity = np.array([100.0, 80.0, 90.0])

    filtered_mz, _, _ = filter_peaks_by_isotopes(
        mz_data,
        intensity,
        np.array([3]),
        mz_tolerance=1e-5,
        units="da",
        use_decreasing_model=use_decreasing_model,
    )

    assert filtered_mz.size == expected_length


@pytest.mark.parametrize(
    ("min_isotope_count", "expected_length"),
    [(2, 1), (3, 2)],
)
def test_requires_minimum_isotope_count(
    min_isotope_count: int,
    expected_length: int,
) -> None:
    mz_data = np.array([100.0, 100.0 + DELTA])
    intensity = np.array([100.0, 80.0])

    filtered_mz, _, _ = filter_peaks_by_isotopes(
        mz_data,
        intensity,
        np.array([2]),
        mz_tolerance=1e-5,
        units="da",
        min_isotope_count=min_isotope_count,
    )

    assert filtered_mz.size == expected_length


def test_maximum_isotope_count_limits_cluster_extension() -> None:
    mz_data = np.array([100.0, 100.0 + DELTA, 100.0 + 2 * DELTA, 100.0 + 3 * DELTA])
    intensity = np.array([100.0, 80.0, 60.0, 40.0])

    filtered_mz, _, _ = filter_peaks_by_isotopes(
        mz_data,
        intensity,
        np.array([4]),
        mz_tolerance=1e-5,
        units="da",
        max_isotope_count=3,
    )

    np.testing.assert_allclose(filtered_mz, [100.0, 100.0 + 3 * DELTA])


# ---------------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------------


def test_rejects_shared_mz_axis() -> None:
    with pytest.raises(
        ValueError,
        match=r"requires processed \+ centroided m/z data",
    ):
        filter_peaks_by_isotopes(
            np.array([100.0, 100.0 + DELTA, 100.0 + 2 * DELTA]),
            np.ones(6),
            np.array([3, 3]),
        )


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"mz_tolerance": 0.0}, "mz_tolerance must be a finite positive number"),
        ({"mz_tolerance": np.nan}, "mz_tolerance must be a finite positive number"),
        (
            {"mz_tolerance": 101.0},
            "mz_tolerance must not exceed 100 ppm or 0.1 Da",
        ),
        (
            {"mz_tolerance": 0.11, "units": "da"},
            "mz_tolerance must not exceed 100 ppm or 0.1 Da",
        ),
        ({"units": "invalid"}, "units must be either 'da' or 'ppm'"),
        ({"min_charge": 0}, "Charges must satisfy"),
        ({"min_charge": 2, "max_charge": 1}, "Charges must satisfy"),
        ({"min_isotope_count": 1}, "Isotope counts must satisfy"),
        (
            {"min_isotope_count": 4, "max_isotope_count": 3},
            "Isotope counts must satisfy",
        ),
        ({"start_intensity_check": -1}, "start_intensity_check must be non-negative"),
    ],
)
def test_rejects_invalid_algorithm_parameters(
    kwargs: dict[str, object],
    message: str,
) -> None:
    with pytest.raises(ValueError, match=message):
        filter_peaks_by_isotopes(
            np.array([100.0, 100.0 + DELTA, 100.0 + 2 * DELTA]),
            np.array([100.0, 80.0, 60.0]),
            np.array([3]),
            **kwargs,  # type: ignore[arg-type]
        )


def test_rejects_lengths_that_do_not_cover_intensity_data() -> None:
    with pytest.raises(ValueError, match=r"sum\(lengths\) must equal intensity.size"):
        filter_peaks_by_isotopes(
            np.array([100.0, 100.0 + DELTA, 100.0 + 2 * DELTA]),
            np.array([100.0, 80.0, 60.0]),
            np.array([2]),
        )
