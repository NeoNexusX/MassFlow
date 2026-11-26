"""Focused contracts for the source-backed MSI adduct filter."""

import csv
from dataclasses import replace
from types import SimpleNamespace
from typing import cast

import matplotlib
from matplotlib.collections import PathCollection
import matplotlib.pyplot as plt
import numpy as np
import pytest

from massflow.data_manager import MSDataManagerImzML
import massflow.preprocess.helper.adduct_filter_helper as adduct_helper
from massflow.preprocess.helper.adduct_filter_helper import (
    CAMERA_COMMIT,
    AdductGroup,
    AdductHypothesis,
    AdductRule,
    PRIMARY_ADDUCT_RULES_NEGATIVE,
    PRIMARY_ADDUCT_RULES_POSITIVE,
    _build_adduct_annotations,
    _build_adduct_keep_mask,
    _build_hypotheses,
    _candidate_groups_from_clusters,
    _cluster_hypotheses,
    _collect_pair_correlations,
    _confirm_groups_from_correlations,
    _neutral_mass,
    _resolve_groups_by_ips,
    _resolve_rules,
    _select_summary_annotation_indices,
    apply_adduct_annotations,
    compute_adduct_annotations,
    write_adduct_filter_report,
)
from massflow.tools.plot_adduct_filter import (
    _plot_correlation_distribution,
    _plot_removal_summary,
    _summary_title,
    plot_adduct_filter_result,
)

matplotlib.use("Agg", force=True)


class _FakeContinuousDataManager:  # pylint: disable=too-few-public-methods
    """Small continuous-data manager exposing the real flat-generator contract."""

    def __init__(
        self,
        mz_axis: np.ndarray,
        batches: list[tuple[np.ndarray, np.ndarray]],
    ) -> None:
        self.ms = SimpleNamespace(
            meta=SimpleNamespace(continuous=True),
            shared_mz_list=np.asarray(mz_axis),
        )
        self.batches = batches
        self.requested_batch_sizes: list[int] = []
        self.target_locs: tuple[list[int], list[int]] | None = None

    def flat_generator(
        self,
        *,
        batch_size: int = 256,
        include_mz: bool = True,
        **_: object,
    ):
        self.requested_batch_sizes.append(batch_size)
        for intensity_flat, lengths in self.batches:
            coordinates = np.zeros((lengths.size, 3), dtype=np.int32)
            yield (
                self.ms.shared_mz_list if include_mz else None,
                intensity_flat,
                lengths,
                coordinates,
            )


def _flat_batch(rows: list[list[float]]) -> tuple[np.ndarray, np.ndarray]:
    lengths = np.asarray([len(row) for row in rows], dtype=np.uint32)
    intensity_flat = np.concatenate([np.asarray(row, dtype=np.float64) for row in rows])
    return intensity_flat, lengths


def _group(
    peaks: tuple[int, ...],
    rules: tuple[int, ...],
    *,
    representative: int,
    mass: float,
    correlations: tuple[float, ...] | None = None,
) -> AdductGroup:
    if correlations is None:
        correlations = tuple(1.0 for _ in peaks)
    return AdductGroup(
        neutral_mass=mass,
        peak_indices=peaks,
        rule_indices=rules,
        neutral_masses=tuple(mass for _ in peaks),
        representative_peak=representative,
        support_peak_indices=tuple(representative for _ in peaks),
        spatial_correlations=correlations,
    )


# ---------------------------------------------------------------------------
# CAMERA 1.68.0 source parity
# ---------------------------------------------------------------------------


def test_primary_rule_rows_match_camera_1_68() -> None:
    assert CAMERA_COMMIT == "2bdda3fecf94e482e86beebb7f981ffeaa9c8e9f"
    assert [
        (r.name, r.nmol, r.charge, r.massdiff, r.oidscore, r.quasi, r.ips)
        for r in PRIMARY_ADDUCT_RULES_POSITIVE
    ] == [
        ("[M+H]+", 1, 1, 1.007276, 1, True, 1.0),
        ("[M+Na]+", 1, 1, 22.989218, 8, True, 1.0),
        ("[M+K]+", 1, 1, 38.963158, 10, True, 1.0),
        ("[M+NH4]+", 1, 1, 18.033823, 16, True, 1.0),
    ]
    assert [
        (r.name, r.nmol, r.charge, r.massdiff, r.oidscore, r.quasi, r.ips)
        for r in PRIMARY_ADDUCT_RULES_NEGATIVE
    ] == [
        ("[M-H]-", 1, -1, -1.007276, 1, True, 1.0),
        ("[M-2H+Na]-", 1, -1, 20.974666, 4, False, 0.5),
        ("[M-2H+K]-", 1, -1, 36.948606, 6, False, 0.5),
        ("[M+Cl]-", 1, -1, 34.969402, 8, True, 1.0),
    ]


def test_neutral_mass_formula_round_trips_all_primary_rules() -> None:
    neutral_mass = 100.0
    for rule in PRIMARY_ADDUCT_RULES_POSITIVE + PRIMARY_ADDUCT_RULES_NEGATIVE:
        observed_mz = (neutral_mass * rule.nmol + rule.massdiff) / abs(rule.charge)
        assert np.isclose(_neutral_mass(observed_mz, rule), neutral_mass)


def test_complete_linkage_rejects_old_moving_mean_chain() -> None:
    hypotheses = [
        AdductHypothesis(0, 0, 100.000),
        AdductHypothesis(1, 1, 100.014),
        AdductHypothesis(2, 2, 100.021),
    ]

    clusters = _cluster_hypotheses(hypotheses, ppm=1e-12, mzabs=0.015)

    # R: cutree(hclust(dist(c(100,100.014,100.021))), h=.015) -> 1,2,2.
    assert [[item.peak_index for item in cluster] for cluster in clusters] == [[1, 2]]


def test_camera_adjacent_integer_bucket_boundary_is_joined() -> None:
    hypotheses = [
        AdductHypothesis(0, 0, 99.995),
        AdductHypothesis(1, 1, 100.005),
    ]

    clusters = _cluster_hypotheses(hypotheses, ppm=1e-12, mzabs=0.015)

    assert len(clusters) == 1
    assert {item.peak_index for item in clusters[0]} == {0, 1}


def test_general_charge_or_molecule_rules_are_rejected_until_camera_checks_exist() -> (
    None
):
    unsupported = (AdductRule("[2M+H]+", 2, 1, 1.007276, 1, True, 1.0),)

    with pytest.raises(ValueError, match="nmol=1"):
        _resolve_rules("positive", unsupported)


# ---------------------------------------------------------------------------
# Candidate selection, spatial network, and CAMERA reduction
# ---------------------------------------------------------------------------


def test_candidate_group_uses_maximum_summed_intensity_representative() -> None:
    mz_axis = np.array(
        [101.007276, 118.033823, 122.989218, 138.963158],
        dtype=np.float64,
    )
    summed = np.array([100.0, 20.0, 400.0, 10.0], dtype=np.float64)
    clusters = _cluster_hypotheses(
        _build_hypotheses(mz_axis, PRIMARY_ADDUCT_RULES_POSITIVE),
        ppm=5.0,
        mzabs=0.015,
    )

    groups = _candidate_groups_from_clusters(
        clusters,
        PRIMARY_ADDUCT_RULES_POSITIVE,
        summed,
        require_quasi=True,
    )

    assert len(groups) == 1
    assert groups[0].representative_peak == 2  # [M+Na]+, not hard-coded [M+H]+


def test_negative_optional_only_group_fails_quasi_requirement() -> None:
    mz_axis = np.array([120.974666, 136.948606], dtype=np.float64)
    summed = np.array([40.0, 20.0], dtype=np.float64)
    clusters = _cluster_hypotheses(
        _build_hypotheses(mz_axis, PRIMARY_ADDUCT_RULES_NEGATIVE),
        ppm=5.0,
        mzabs=0.015,
    )

    strict = _candidate_groups_from_clusters(
        clusters,
        PRIMARY_ADDUCT_RULES_NEGATIVE,
        summed,
        require_quasi=True,
    )
    relaxed = _candidate_groups_from_clusters(
        clusters,
        PRIMARY_ADDUCT_RULES_NEGATIVE,
        summed,
        require_quasi=False,
    )

    assert not strict
    assert len(relaxed) == 1


def test_spatial_graph_splits_two_molecules_in_one_mass_cluster() -> None:
    candidate = AdductGroup(
        neutral_mass=500.002,
        peak_indices=(0, 1, 2, 3),
        rule_indices=(0, 1, 0, 1),
        neutral_masses=(500.000, 500.000, 500.004, 500.004),
        representative_peak=0,
        support_peak_indices=(0, 1, 2, 3),
        spatial_correlations=(np.nan, np.nan, np.nan, np.nan),
    )
    correlations = {
        (0, 1): 0.99,
        (0, 3): -0.90,
        (1, 2): -0.90,
        (2, 3): 0.98,
    }

    confirmed, evidence = _confirm_groups_from_correlations(
        [candidate],
        correlations,
        np.array([100.0, 40.0, 80.0, 30.0]),
        PRIMARY_ADDUCT_RULES_POSITIVE,
        require_quasi=True,
        min_spatial_correlation=0.7,
    )

    assert [group.peak_indices for group in confirmed] == [(0, 1), (2, 3)]
    assert [group.representative_peak for group in confirmed] == [0, 2]
    assert sum(item.passes_threshold for item in evidence) == 2


def test_correlation_threshold_is_inclusive() -> None:
    candidate = _group((0, 1), (0, 1), representative=0, mass=100.0)

    confirmed, _ = _confirm_groups_from_correlations(
        [candidate],
        {(0, 1): 0.7},
        np.array([10.0, 5.0]),
        PRIMARY_ADDUCT_RULES_POSITIVE,
        require_quasi=True,
        min_spatial_correlation=0.7,
    )
    rejected, _ = _confirm_groups_from_correlations(
        [candidate],
        {(0, 1): np.nextafter(0.7, 0.0)},
        np.array([10.0, 5.0]),
        PRIMARY_ADDUCT_RULES_POSITIVE,
        require_quasi=True,
        min_spatial_correlation=0.7,
    )

    assert len(confirmed) == 1
    assert not rejected


def test_ips_discards_lower_scoring_overlapping_mass_group() -> None:
    lower = _group((0, 1), (0, 1), representative=0, mass=100.0)
    higher = _group((0, 2, 3), (0, 1, 2), representative=0, mass=200.0)

    resolved = _resolve_groups_by_ips(
        [lower, higher],
        PRIMARY_ADDUCT_RULES_POSITIVE,
    )

    assert resolved == [higher]


# ---------------------------------------------------------------------------
# Streaming Pearson and preview/apply contract
# ---------------------------------------------------------------------------


def test_streaming_pair_correlations_match_concatenated_numpy_reference() -> None:
    mz_axis = np.arange(4, dtype=np.float64)
    batches = [
        _flat_batch([[1.0, 2.0, 3.0, 4.0], [2.0, 4.0, 6.0]]),
        _flat_batch([[3.0, 6.0, 9.0, 12.0], [4.0, 8.0]]),
    ]
    manager = _FakeContinuousDataManager(mz_axis, batches)
    pairs = [(0, 1), (0, 3), (2, 3)]

    actual = _collect_pair_correlations(
        cast(MSDataManagerImzML, manager),
        pairs,
        batch_size=2,
    )

    padded = np.array(
        [
            [1.0, 2.0, 3.0, 4.0],
            [2.0, 4.0, 6.0, 0.0],
            [3.0, 6.0, 9.0, 12.0],
            [4.0, 8.0, 0.0, 0.0],
        ]
    )
    for pair in pairs:
        expected = float(np.corrcoef(padded[:, pair[0]], padded[:, pair[1]])[0, 1])
        assert np.isclose(actual[pair], expected)
    assert manager.requested_batch_sizes == [2]


def test_public_preview_does_not_remove_when_no_rule_pair_matches() -> None:
    manager = _FakeContinuousDataManager(
        np.array([101.007276, 250.0]),
        [_flat_batch([[10.0, 4.0], [20.0, 8.0], [30.0, 12.0]])],
    )

    result = compute_adduct_annotations(manager, batch_size=3)  # type: ignore[arg-type]

    assert not result.groups
    assert not result.annotations
    assert not result.pair_evidence


def test_public_preview_rejects_rule_match_with_low_spatial_correlation() -> None:
    manager = _FakeContinuousDataManager(
        np.array([101.007276, 122.989218]),
        [
            _flat_batch(
                [
                    [10.0, 40.0],
                    [20.0, 30.0],
                    [30.0, 20.0],
                    [40.0, 10.0],
                ]
            )
        ],
    )

    result = compute_adduct_annotations(manager, batch_size=4)  # type: ignore[arg-type]

    assert not result.groups
    assert not result.annotations
    assert len(result.pair_evidence) == 1
    assert not result.pair_evidence[0].passes_threshold


def test_public_preview_keeps_two_spatial_components_from_one_mass_cluster() -> None:
    manager = _FakeContinuousDataManager(
        np.array([501.007276, 501.011276, 522.989218, 522.993218]),
        [
            _flat_batch(
                [
                    [100.0, 0.0, 40.0, 0.0],
                    [80.0, 0.0, 30.0, 0.0],
                    [0.0, 70.0, 0.0, 25.0],
                    [0.0, 90.0, 0.0, 35.0],
                ]
            )
        ],
    )

    result = compute_adduct_annotations(manager, batch_size=4)  # type: ignore[arg-type]

    assert [group.peak_indices for group in result.groups] == [(0, 2), (1, 3)]
    assert {annotation.removed_peak for annotation in result.annotations} == {2, 3}


@pytest.fixture(name="synthetic_manager")
def _synthetic_manager() -> _FakeContinuousDataManager:
    mz_axis = np.array(
        [101.007276, 118.033823, 122.989218, 138.963158, 250.0],
        dtype=np.float64,
    )
    batches = [
        _flat_batch(
            [
                [10.0, 5.0, 20.0, 2.0, 1.0],
                [20.0, 10.0, 40.0, 4.0, 2.0],
            ]
        ),
        _flat_batch(
            [
                [30.0, 15.0, 60.0, 6.0, 4.0],
                [40.0, 20.0, 80.0, 8.0, 3.0],
            ]
        ),
    ]
    return _FakeContinuousDataManager(mz_axis, batches)


def test_compute_annotations_returns_immutable_auditable_preview(
    synthetic_manager: _FakeContinuousDataManager,
) -> None:
    result = compute_adduct_annotations(
        synthetic_manager,  # type: ignore[arg-type]
        batch_size=2,
    )

    assert len(result.groups) == 1
    assert result.groups[0].representative_peak == 2
    assert {annotation.removed_peak for annotation in result.annotations} == {0, 1, 3}
    assert all(
        annotation.spatial_correlation >= 0.7 for annotation in result.annotations
    )
    assert result.config.min_spatial_correlation == 0.7
    assert result.config.ppm == 5.0
    assert len(result.pair_evidence) == 6
    with pytest.raises(ValueError):
        result.mz_axis[0] = 0.0


def test_apply_reuses_preview_mask_and_writes_matching_report(
    synthetic_manager: _FakeContinuousDataManager,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path,
) -> None:
    result = compute_adduct_annotations(
        synthetic_manager,  # type: ignore[arg-type]
        batch_size=2,
    )
    captured: dict[str, np.ndarray] = {}
    sentinel = object()

    def fake_filter(_manager, keep_mask, **_kwargs):
        captured["keep_mask"] = keep_mask.copy()
        return sentinel

    monkeypatch.setattr(adduct_helper, "filter_data_by_mask", fake_filter)
    report_base_path = tmp_path / "adduct-report"
    report_path = report_base_path.with_suffix(".csv")

    returned = apply_adduct_annotations(
        synthetic_manager,  # type: ignore[arg-type]
        result,
        report_path=report_base_path,
    )

    assert returned is sentinel
    np.testing.assert_array_equal(
        captured["keep_mask"],
        np.array([False, False, True, False, True]),
    )
    with report_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        assert reader.fieldnames == [
            "representative_mz", "adduct_mz", "representative_adduct", "adduct"
        ]
        rows = list(reader)
    selected = _select_summary_annotation_indices(result.annotations)
    annotations = [
        annotation for index, annotation in enumerate(result.annotations)
        if index in selected
    ]
    assert len(rows) == len({item.removed_peak for item in result.annotations})
    for row, annotation in zip(rows, annotations, strict=True):
        assert float(row["representative_mz"]) == annotation.representative_mz
        assert float(row["adduct_mz"]) == annotation.removed_mz
        assert row["representative_adduct"] == annotation.representative_rule
        assert row["adduct"] == annotation.removed_rule


def test_adduct_report_defaults_to_no_removal(
    synthetic_manager: _FakeContinuousDataManager,
    tmp_path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def unexpected_removal(*_args, **_kwargs):
        pytest.fail("reporting must not apply a peak mask")

    monkeypatch.setattr(adduct_helper, "filter_data_by_mask", unexpected_removal)
    report_path = tmp_path / "preview.csv"
    returned = adduct_helper.filter_peaks_by_adducts(
        synthetic_manager,  # type: ignore[arg-type]
        report_path=report_path,
    )
    assert returned is synthetic_manager
    with report_path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 3


def test_report_contains_one_row_per_peak_and_handles_no_matches(
    synthetic_manager: _FakeContinuousDataManager,
    tmp_path,
) -> None:
    result = compute_adduct_annotations(synthetic_manager)  # type: ignore[arg-type]
    duplicated = replace(
        result, groups=result.groups * 2, annotations=result.annotations * 2,
    )
    report_path = write_adduct_filter_report(duplicated, tmp_path / "duplicates.csv")
    with report_path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == len({item.removed_peak for item in result.annotations})
    assert len({row["adduct_mz"] for row in rows}) == len(rows)

    empty = replace(result, groups=(), annotations=(), pair_evidence=())
    empty_path = write_adduct_filter_report(empty, tmp_path / "empty.csv")
    assert empty_path.read_text(encoding="utf-8").splitlines() == [
        "representative_mz,adduct_mz,representative_adduct,adduct"
    ]


def test_report_rejects_non_csv_extension(
    synthetic_manager: _FakeContinuousDataManager,
    tmp_path,
) -> None:
    result = compute_adduct_annotations(
        synthetic_manager,  # type: ignore[arg-type]
        batch_size=2,
    )

    with pytest.raises(ValueError, match=".csv extension"):
        write_adduct_filter_report(result, tmp_path / "adduct-report.txt")


def test_apply_rejects_result_with_threshold_inconsistent_annotation(
    synthetic_manager: _FakeContinuousDataManager,
) -> None:
    result = compute_adduct_annotations(
        synthetic_manager,  # type: ignore[arg-type]
        batch_size=2,
    )
    inconsistent = replace(
        result,
        annotations=(
            replace(result.annotations[0], correlation_threshold=0.9),
            *result.annotations[1:],
        ),
    )

    with pytest.raises(ValueError, match="threshold does not match"):
        apply_adduct_annotations(
            synthetic_manager,  # type: ignore[arg-type]
            inconsistent,
        )


def test_apply_rejects_different_manager_with_identical_mz_axis(
    synthetic_manager: _FakeContinuousDataManager,
) -> None:
    result = compute_adduct_annotations(
        synthetic_manager,  # type: ignore[arg-type]
        batch_size=2,
    )
    other_manager = _FakeContinuousDataManager(
        synthetic_manager.ms.shared_mz_list.copy(),
        list(synthetic_manager.batches),
    )

    with pytest.raises(ValueError, match="different data manager"):
        apply_adduct_annotations(
            other_manager,  # type: ignore[arg-type]
            result,
        )


def test_apply_rejects_target_locs_changed_after_preview(
    synthetic_manager: _FakeContinuousDataManager,
) -> None:
    synthetic_manager.target_locs = ([0, 0], [10, 10])
    result = compute_adduct_annotations(
        synthetic_manager,  # type: ignore[arg-type]
        batch_size=2,
    )
    synthetic_manager.target_locs = ([0, 0], [5, 5])

    with pytest.raises(ValueError, match="target_locs"):
        apply_adduct_annotations(
            synthetic_manager,  # type: ignore[arg-type]
            result,
        )


def test_annotations_and_keep_mask_remove_the_same_unique_peaks() -> None:
    groups = [
        _group((0, 1), (0, 1), representative=1, mass=100.0),
        _group((1, 2, 3), (1, 2, 3), representative=3, mass=200.0),
    ]
    mz_axis = np.array([100.0, 200.0, 300.0, 400.0])

    annotations = _build_adduct_annotations(
        groups,
        mz_axis,
        PRIMARY_ADDUCT_RULES_POSITIVE,
        min_spatial_correlation=0.7,
    )
    keep = _build_adduct_keep_mask(mz_axis.size, groups)

    assert {annotation.removed_peak for annotation in annotations} == set(
        np.flatnonzero(~keep)
    )
    assert keep[1]  # representative in the first group remains protected


def test_summary_evidence_selects_strongest_row_per_removed_peak(
    synthetic_manager: _FakeContinuousDataManager,
) -> None:
    result = compute_adduct_annotations(
        synthetic_manager,  # type: ignore[arg-type]
        batch_size=2,
    )
    stronger = result.annotations[0]
    weaker = replace(stronger, spatial_correlation=0.7)
    annotations = (weaker, stronger, *result.annotations[1:])

    selected = _select_summary_annotation_indices(annotations)

    assert 0 not in selected
    assert 1 in selected
    assert len(selected) == len({item.removed_peak for item in annotations})


# ---------------------------------------------------------------------------
# Decision-preview rendering
# ---------------------------------------------------------------------------


def test_plot_uses_recorded_threshold_and_unique_removed_count(
    synthetic_manager: _FakeContinuousDataManager,
) -> None:
    result = compute_adduct_annotations(
        synthetic_manager,  # type: ignore[arg-type]
        batch_size=2,
    )
    duplicated = replace(
        result,
        annotations=result.annotations + (result.annotations[0],),
    )

    title = _summary_title(duplicated)
    assert "5 peaks" in title
    assert "3 unique peaks removed" in title
    assert "4 evidence rows" in title

    figure, axis = plt.subplots()
    _plot_correlation_distribution(axis, result)
    np.testing.assert_allclose(
        np.asarray(axis.lines[0].get_xdata(), dtype=np.float64),
        result.config.min_spatial_correlation,
    )
    plt.close(figure)


def test_removal_summary_plots_one_point_per_unique_removed_peak(
    synthetic_manager: _FakeContinuousDataManager,
) -> None:
    result = compute_adduct_annotations(
        synthetic_manager,  # type: ignore[arg-type]
        batch_size=2,
    )
    figure, axis = plt.subplots()

    _plot_removal_summary(axis, result)

    plotted_points = sum(
        cast(PathCollection, collection).get_offsets().shape[0] # type: ignore
        for collection in axis.collections
    )
    assert plotted_points == len(
        {annotation.removed_peak for annotation in result.annotations}
    )
    plt.close(figure)


def test_plot_result_writes_png_and_closes_figure(
    synthetic_manager: _FakeContinuousDataManager,
    tmp_path,
) -> None:
    result = compute_adduct_annotations(
        synthetic_manager,  # type: ignore[arg-type]
        batch_size=2,
    )
    output_path = tmp_path / "adduct-preview.png"
    report_base_path = tmp_path / "adduct-preview-report"
    figures_before = tuple(plt.get_fignums())

    returned = plot_adduct_filter_result(
        result,
        show=False,
        output_path=output_path,
        report_path=report_base_path,
    )

    assert returned is result
    assert output_path.stat().st_size > 0
    with report_base_path.with_suffix(".csv").open(
        newline="", encoding="utf-8"
    ) as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == len({item.removed_peak for item in result.annotations})
    assert tuple(plt.get_fignums()) == figures_before


def test_empty_plot_result_is_renderable(
    synthetic_manager: _FakeContinuousDataManager,
    tmp_path,
) -> None:
    result = compute_adduct_annotations(
        synthetic_manager,  # type: ignore[arg-type]
        batch_size=2,
    )
    empty = replace(result, annotations=(), groups=(), pair_evidence=())
    output_path = tmp_path / "empty-adduct-preview.png"

    plot_adduct_filter_result(empty, show=False, output_path=output_path)

    assert output_path.stat().st_size > 0


def test_plot_rejects_annotation_rows_not_reconstructible_from_groups(
    synthetic_manager: _FakeContinuousDataManager,
    tmp_path,
) -> None:
    result = compute_adduct_annotations(
        synthetic_manager,  # type: ignore[arg-type]
        batch_size=2,
    )
    duplicated = replace(
        result,
        annotations=result.annotations + (result.annotations[0],),
    )

    with pytest.raises(ValueError, match="rows do not match"):
        plot_adduct_filter_result(
            duplicated,
            show=False,
            output_path=tmp_path / "must-not-render.png",
        )
