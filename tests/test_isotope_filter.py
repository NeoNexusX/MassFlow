"""Algorithm-level tests for isotope peak filtering."""

import csv

import numpy as np
import pytest

from massflow.data_manager import MSDataManagerImzML
from massflow.preprocess.flat_pre_fun import FlatPreprocess
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
        remove=True,
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
        remove=True,
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
        remove=True,
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
        remove=True,
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
        remove=True,
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
        remove=True,
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
        remove=True,
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
        remove=True,
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
        remove=True,
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
        remove=True,
        mz_tolerance=1e-5,
        units="da",
        max_isotope_count=3,
    )

    np.testing.assert_allclose(filtered_mz, [100.0, 100.0 + 3 * DELTA])


# ---------------------------------------------------------------------------
# Reporting without removal
# ---------------------------------------------------------------------------


def test_default_reports_each_spectrum_without_removing_peaks(tmp_path) -> None:
    mz_data = np.array(
        [100.0, 100.0 + DELTA, 100.0 + 2 * DELTA,
         500.0, 500.0 + DELTA / 2, 500.0 + DELTA]
    )
    intensity = np.array([100.0, 80.0, 60.0, 90.0, 70.0, 50.0])
    lengths = np.array([3, 3], dtype=np.uint32)
    report_path = tmp_path / "isotope-report"

    mz, values, output_lengths = filter_peaks_by_isotopes(
        mz_data=mz_data, intensity=intensity, lengths=lengths,
        max_charge=2,
        report_path=report_path,
    )

    assert mz is mz_data
    assert values is intensity
    np.testing.assert_array_equal(output_lengths, lengths)
    with report_path.with_suffix(".csv").open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        assert reader.fieldnames == [
            "spectrum_index", "representative_mz", "isotope_mz", "charge"
        ]
        rows = list(reader)
    assert [int(row["spectrum_index"]) for row in rows] == [0, 0, 1, 1]
    assert [int(row["charge"]) for row in rows] == [1, 1, 2, 2]
    np.testing.assert_allclose(
        [float(row["representative_mz"]) for row in rows],
        [100.0, 100.0, 500.0, 500.0],
    )
    np.testing.assert_allclose(
        [float(row["isotope_mz"]) for row in rows], mz_data[[1, 2, 4, 5]]
    )

    removal_report = tmp_path / "isotope-removal.csv"
    filtered_mz, _, filtered_lengths = filter_peaks_by_isotopes(
        mz_data, intensity, lengths,
        max_charge=2,
        remove=True,
        report_path=removal_report,
    )
    np.testing.assert_allclose(filtered_mz, [100.0, 500.0])
    np.testing.assert_array_equal(filtered_lengths, [1, 1])
    assert removal_report.read_bytes() == report_path.with_suffix(".csv").read_bytes()


def test_report_without_matches_contains_only_header(tmp_path) -> None:
    report_path = tmp_path / "empty.csv"
    filter_peaks_by_isotopes(
        np.array([100.0, 150.0]), np.array([5.0, 2.0]), np.array([2]),
        report_path=report_path,
    )
    assert report_path.read_text().splitlines() == [
        "spectrum_index,representative_mz,isotope_mz,charge"
    ]


def test_flat_task_always_removes_isotope_peaks() -> None:
    mz_data = np.array([100.0, 100.0 + DELTA, 100.0 + 2 * DELTA, 150.0])
    intensity = np.array([100.0, 80.0, 60.0, 5.0])

    result = FlatPreprocess.isotope_filter_flat(
        mz_data, intensity, np.array([4]),
    )

    assert result.mz_data is not None
    np.testing.assert_allclose(result.mz_data, [100.0, 150.0])
    np.testing.assert_allclose(result.intensity, [100.0, 5.0])
    np.testing.assert_array_equal(result.lengths, [2])


# ---------------------------------------------------------------------------
# imzML data-manager input
# ---------------------------------------------------------------------------


@pytest.fixture(name="processed_imzml")
def processed_imzml_fixture(tmp_path):
    with MSDataManagerImzML(filepath=str(tmp_path / "source.imzML")) as source:
        source.ms.meta.processed = True
        source.ms.meta.centroid_spectrum = True
        source.write_flat_batch(
            mz_axis=np.array([
                80.0, 90.0,
                100.0, 100.0 + DELTA, 100.0 + 2 * DELTA,
                500.0, 500.0 + DELTA / 2, 500.0 + DELTA,
            ]),
            intensity_flat=np.array(
                [5.0, 4.0, 100.0, 80.0, 60.0, 90.0, 70.0, 50.0],
                dtype=np.float32,
            ),
            lengths=np.array([2, 3, 3], dtype=np.uint32),
            coordinates=np.array([[1, 1, 1], [2, 1, 1], [3, 1, 1]]),
        )
        source.close_writer()
        source.load_head_data()
        yield source


def test_imzml_report_only_streams_all_batches_without_creating_output(
    processed_imzml, tmp_path,
) -> None:
    report_path = tmp_path / "report"
    output_path = tmp_path / "unused.imzML"
    returned = filter_peaks_by_isotopes(
        data_manager=processed_imzml,
        max_charge=2,
        batch_size=1,
        report_path=report_path,
        output_path=output_path,
    )

    assert returned is processed_imzml
    assert not output_path.exists()
    assert not output_path.with_suffix(".ibd").exists()
    with report_path.with_suffix(".csv").open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        assert reader.fieldnames == [
            "spectrum_index", "representative_mz", "isotope_mz", "charge"
        ]
        rows = list(reader)
    assert [int(row["spectrum_index"]) for row in rows] == [1, 1, 2, 2]
    assert [int(row["charge"]) for row in rows] == [1, 1, 2, 2]
    np.testing.assert_allclose(
        [float(row["representative_mz"]) for row in rows],
        [100.0, 100.0, 500.0, 500.0],
    )
    _, _, lengths, _ = next(processed_imzml.flat_generator(batch_size=3))
    np.testing.assert_array_equal(lengths, [2, 3, 3])


@pytest.mark.parametrize("temporary_output", [False, True])
def test_imzml_removal_returns_readable_filtered_manager(
    processed_imzml, tmp_path, temporary_output: bool,
) -> None:
    output_path = None if temporary_output else tmp_path / "filtered.imzML"
    report_path = tmp_path / "filtered.csv"
    filtered = filter_peaks_by_isotopes(
        data_manager=processed_imzml,
        max_charge=2,
        batch_size=1,
        remove=True,
        output_path=output_path,
        temp_dir=tmp_path,
        report_path=report_path,
    )

    assert isinstance(filtered, MSDataManagerImzML)
    assert filtered is not processed_imzml
    with filtered:
        assert filtered.ms.meta.processed is True
        assert filtered.ms.meta.continuous is not True
        assert filtered.ms.meta.centroid_spectrum is True
        assert filtered.ms.meta.profile_spectrum is not True
        mz, values, lengths, coordinates = next(filtered.flat_generator(batch_size=3))
        np.testing.assert_allclose(mz, [80.0, 90.0, 100.0, 500.0])
        np.testing.assert_allclose(values, [5.0, 4.0, 100.0, 90.0])
        np.testing.assert_array_equal(lengths, [2, 1, 1])
        np.testing.assert_array_equal(coordinates, [[1, 1, 1], [2, 1, 1], [3, 1, 1]])
        assert mz.dtype == np.dtype(processed_imzml.mz_dtype)
        assert values.dtype == np.dtype(processed_imzml.intensity_dtype)
    _, _, source_lengths, _ = next(processed_imzml.flat_generator(batch_size=3))
    np.testing.assert_array_equal(source_lengths, [2, 3, 3])
    with report_path.open(newline="", encoding="utf-8") as handle:
        assert len(list(csv.DictReader(handle))) == 4


def test_imzml_removal_rejects_overwriting_source(processed_imzml) -> None:
    with pytest.raises(ValueError, match="must differ from the input"):
        filter_peaks_by_isotopes(
            data_manager=processed_imzml,
            remove=True,
            output_path=processed_imzml.imzml_filepath,
        )
    _, _, lengths, _ = next(processed_imzml.flat_generator(batch_size=3))
    np.testing.assert_array_equal(lengths, [2, 3, 3])


@pytest.mark.parametrize("metadata_field", ["continuous", "profile_spectrum"])
def test_imzml_input_requires_processed_centroided_spectra(
    processed_imzml, metadata_field: str,
) -> None:
    setattr(processed_imzml.ms.meta, metadata_field, True)
    if metadata_field == "profile_spectrum":
        processed_imzml.ms.meta.centroid_spectrum = None
    with pytest.raises(ValueError, match=r"processed \+ centroided imzML"):
        filter_peaks_by_isotopes(data_manager=processed_imzml)


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
