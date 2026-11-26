"""Source-backed adduct annotation and conservative reduction for MSI data.

The source-backed part of this module follows CAMERA 1.68.0:

* the bundled primary-rule fields and values come from
  ``inst/rules/primary_adducts_{pos,neg}.csv``;
* neutral masses use CAMERA ``massDiffMatrix``;
* hypothesis clustering follows ``createHypothese`` (integer mass buckets,
  adjacent-bucket boundary candidates, and complete-linkage clustering);
* quasi-ion and IPS competition checks follow ``checkQuasimolion`` and
  ``checkIps``; and
* reduction keeps the highest summed-intensity feature, matching
  ``getReducedPeaklist(..., method="maxint", default.adduct.info="maxint")``.

MassFlow/MSI adaptations are deliberately separate:

* the four CAMERA primary CSV rules are a narrow MassFlow preset. CAMERA's
  ``rules=NULL`` dynamically creates a larger rule table;
* because MSI has no LC-MS retention-time pseudospectra, neutral-mass
  candidates are generated on the shared m/z axis and then split into
  connected components of spatially correlated candidate pairs; and
* ``min_spatial_correlation`` is a MassFlow threshold. The CAMERA defaults
  ``ppm=5`` and ``mzabs=0.015`` are copied for traceability, but their
  suitability for a particular MSI instrument must be assessed by the user.

This first implementation intentionally accepts only singly charged,
single-molecule rules. That covers the bundled primary preset and avoids
claiming parity with CAMERA's isotope-charge and OID-causality handling for
more general rule tables. Isotope filtering should be run before this step.
"""

import csv
from collections import defaultdict
from dataclasses import dataclass, field
from itertools import combinations
from pathlib import Path
from typing import Literal

import numpy as np
from scipy.cluster.hierarchy import fclusterdata

from massflow.data_manager import MSDataManagerImzML
from massflow.preprocess.helper.stat_filter_helper import (
    compute_sum_intensity,
    filter_data_by_mask,
)
from massflow.preprocess.numba.adduct_filter_numba import (
    update_pair_statistics_flat_parallel_jit,
)

Polarity = Literal["positive", "negative"]

CAMERA_VERSION = "1.68.0"
CAMERA_COMMIT = "2bdda3fecf94e482e86beebb7f981ffeaa9c8e9f"
CAMERA_ALGORITHM_SOURCE = (
    "CAMERA 1.68.0 createHypothese/checkIps/getReducedPeaklist; "
    "https://code.bioconductor.org/browse/CAMERA/blob/RELEASE_3_23/"
    "R/fct_findAdducts.R; "
    "https://code.bioconductor.org/browse/CAMERA/blob/RELEASE_3_23/"
    "R/xsAnnotate.R "
    f"(commit {CAMERA_COMMIT})"
)
CAMERA_RULE_SOURCE_POSITIVE = (
    "CAMERA 1.68.0 "
    "https://code.bioconductor.org/browse/CAMERA/blob/RELEASE_3_23/"
    "inst/rules/primary_adducts_pos.csv "
    f"(commit {CAMERA_COMMIT})"
)
CAMERA_RULE_SOURCE_NEGATIVE = (
    "CAMERA 1.68.0 "
    "https://code.bioconductor.org/browse/CAMERA/blob/RELEASE_3_23/"
    "inst/rules/primary_adducts_neg.csv "
    f"(commit {CAMERA_COMMIT})"
)
MSI_SPATIAL_ADAPTATION_SOURCE = (
    "MassFlow MSI adaptation: connected components of accepted spatial "
    "candidate pairs, informed by MZmine IonNetworkingTask correlation-group "
    "networking; https://github.com/mzmine/mzmine/blob/"
    "644a28e697ae46d15b51a2d19abc590ac3e2d725/"
    "mzmine-community/src/main/java/io/github/mzmine/modules/dataprocessing/"
    "id_ion_identity_networking/ionidnetworking/IonNetworkingTask.java"
)


# ---------------------------------------------------------------------------
# CAMERA primary rule preset
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class AdductRule:
    """One CAMERA rule row: name, nmol, charge, massdiff, oidscore, quasi, ips."""

    name: str
    nmol: int
    charge: int
    massdiff: float
    oidscore: int
    quasi: bool
    ips: float


PRIMARY_ADDUCT_RULES_POSITIVE: tuple[AdductRule, ...] = (
    AdductRule("[M+H]+", 1, 1, 1.007276, 1, True, 1.0),
    AdductRule("[M+Na]+", 1, 1, 22.989218, 8, True, 1.0),
    AdductRule("[M+K]+", 1, 1, 38.963158, 10, True, 1.0),
    AdductRule("[M+NH4]+", 1, 1, 18.033823, 16, True, 1.0),
)

PRIMARY_ADDUCT_RULES_NEGATIVE: tuple[AdductRule, ...] = (
    AdductRule("[M-H]-", 1, -1, -1.007276, 1, True, 1.0),
    AdductRule("[M-2H+Na]-", 1, -1, 20.974666, 4, False, 0.5),
    AdductRule("[M-2H+K]-", 1, -1, 36.948606, 6, False, 0.5),
    AdductRule("[M+Cl]-", 1, -1, 34.969402, 8, True, 1.0),
)


def _resolve_rules(
    polarity: Polarity,
    rules: tuple[AdductRule, ...] | None,
) -> tuple[AdductRule, ...]:
    if polarity not in {"positive", "negative"}:
        raise ValueError("polarity must be either 'positive' or 'negative'.")
    selected = (
        (
            PRIMARY_ADDUCT_RULES_POSITIVE
            if polarity == "positive"
            else PRIMARY_ADDUCT_RULES_NEGATIVE
        )
        if rules is None
        else tuple(rules)
    )
    if len({rule.name for rule in selected}) != len(selected):
        raise ValueError("adduct rule names must be unique.")

    expected_charge_sign = 1 if polarity == "positive" else -1
    for rule in selected:
        if rule.nmol != 1 or abs(rule.charge) != 1:
            raise ValueError(
                "This adduct filter currently supports only nmol=1 and "
                "abs(charge)=1 rules; general CAMERA isotope/OID checks are "
                "not implemented."
            )
        if rule.charge != expected_charge_sign:
            raise ValueError(f"rule {rule.name!r} does not match {polarity} polarity.")
        if not np.isfinite(rule.massdiff) or not np.isfinite(rule.ips):
            raise ValueError(f"rule {rule.name!r} has a non-finite massdiff or ips.")
    return selected


# CAMERA neutral-mass hypotheses
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class AdductHypothesis:
    peak_index: int
    rule_index: int
    neutral_mass: float


def _neutral_mass(mz: float, rule: AdductRule) -> float:
    """Return CAMERA/MZmine-equivalent neutral mass for one rule and m/z."""
    return (abs(rule.charge) * mz - rule.massdiff) / rule.nmol


def _build_hypotheses(
    mz_axis: np.ndarray,
    rules: tuple[AdductRule, ...],
) -> list[AdductHypothesis]:
    """Build all feature-rule hypotheses with vectorized neutral-mass arithmetic."""
    mz_values = np.asarray(mz_axis, dtype=np.float64)
    charges = np.asarray([abs(rule.charge) for rule in rules], dtype=np.float64)
    massdiffs = np.asarray([rule.massdiff for rule in rules], dtype=np.float64)
    nmols = np.asarray([rule.nmol for rule in rules], dtype=np.float64)
    neutral = (mz_values[:, None] * charges - massdiffs) / nmols

    neutral_values = neutral.ravel()
    order = np.argsort(neutral_values, kind="stable")
    hypotheses: list[AdductHypothesis] = []
    for index in order:
        if not np.isfinite(neutral_values[index]):
            continue
        peak_index, rule_index = divmod(int(index), len(rules))
        hypotheses.append(
            AdductHypothesis(
                peak_index=peak_index,
                rule_index=rule_index,
                neutral_mass=float(neutral_values[index]),
            )
        )
    return hypotheses


def _neutral_mass_tolerance(
    neutral_mass: float,
    *,
    ppm: float,
    mzabs: float,
) -> float:
    """CAMERA ``createHypothese`` cut height: 2*ppm*mass + mzabs."""
    return 2.0 * ppm * 1e-6 * neutral_mass + mzabs


def _cluster_hypotheses(
    hypotheses: list[AdductHypothesis],
    *,
    ppm: float,
    mzabs: float,
) -> list[list[AdductHypothesis]]:
    """Port CAMERA ``createHypothese`` integer-bucket/complete-linkage grouping."""
    buckets: dict[int, list[AdductHypothesis]] = defaultdict(list)
    for hypothesis in hypotheses:
        integer_mass = int(np.trunc(hypothesis.neutral_mass))
        if integer_mass > 1:
            buckets[integer_mass].append(hypothesis)

    clusters: list[list[AdductHypothesis]] = []
    seen: set[tuple[tuple[int, int], ...]] = set()
    for integer_mass in sorted(buckets):
        candidates = list(buckets[integer_mass])
        boundary_tolerance = _neutral_mass_tolerance(
            float(integer_mass),
            ppm=ppm,
            mzabs=mzabs,
        )
        previous = buckets.get(integer_mass - 1)
        if (
            previous
            and candidates
            and min(item.neutral_mass for item in candidates)
            < integer_mass + boundary_tolerance
        ):
            candidates.extend(
                item
                for item in previous
                if item.neutral_mass > integer_mass - boundary_tolerance
            )
        if len(candidates) < 2:
            continue

        masses = np.asarray(
            [item.neutral_mass for item in candidates],
            dtype=np.float64,
        )
        tolerance = _neutral_mass_tolerance(
            float(np.mean(masses)),
            ppm=ppm,
            mzabs=mzabs,
        )
        labels = fclusterdata(
            masses[:, None],
            t=tolerance,
            criterion="distance",
            metric="euclidean",
            method="complete",
        )
        for label in np.unique(labels):
            member_indices = np.flatnonzero(labels == label)
            if member_indices.size < 2:
                continue
            cluster = [candidates[int(index)] for index in member_indices]
            if len({item.peak_index for item in cluster}) < 2:
                continue
            key = tuple(sorted((item.peak_index, item.rule_index) for item in cluster))
            if key in seen:
                continue
            seen.add(key)
            clusters.append(cluster)
    return clusters


# ---------------------------------------------------------------------------
# Candidate mass groups and spatial evidence
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class AdductGroup:
    neutral_mass: float
    peak_indices: tuple[int, ...]
    rule_indices: tuple[int, ...]
    neutral_masses: tuple[float, ...]
    representative_peak: int
    support_peak_indices: tuple[int, ...]
    spatial_correlations: tuple[float, ...]


@dataclass(frozen=True, slots=True)
class AdductPairEvidence:
    neutral_mass: float
    left_peak: int
    right_peak: int
    left_rule: int
    right_rule: int
    spatial_correlation: float
    passes_threshold: bool


@dataclass(frozen=True, slots=True)
class AdductAnnotation:
    neutral_mass: float
    neutral_mass_span: float
    representative_peak: int
    removed_peak: int
    representative_mz: float
    removed_mz: float
    representative_rule: str
    removed_rule: str
    support_peak: int
    support_mz: float
    support_rule: str
    spatial_correlation: float
    correlation_threshold: float


@dataclass(frozen=True, slots=True)
class AdductFilterConfig:
    polarity: Polarity
    ppm: float
    mzabs: float
    require_quasi: bool
    min_spatial_correlation: float
    batch_size: int


@dataclass(frozen=True, slots=True)
class AdductAnnotationResult:
    """Immutable preview tied to the exact data-manager instance it inspected."""

    mz_axis: np.ndarray
    annotations: tuple[AdductAnnotation, ...]
    groups: tuple[AdductGroup, ...]
    pair_evidence: tuple[AdductPairEvidence, ...]
    rules: tuple[AdductRule, ...]
    rule_names: tuple[str, ...]
    config: AdductFilterConfig
    source_data_manager: object = field(repr=False, compare=False)
    source_manager_state: tuple[object, ...] = field(repr=False, compare=False)


def _choose_representative_peak(
    cluster: list[AdductHypothesis],
    rules: tuple[AdductRule, ...],
    summed_intensity: np.ndarray,
) -> int:
    """Choose CAMERA ``getReducedPeaklist(maxint)``-style representative."""
    del rules
    peak_indices = sorted({item.peak_index for item in cluster})
    return min(
        peak_indices,
        key=lambda peak_index: (-float(summed_intensity[peak_index]), peak_index),
    )


def _choose_rule_per_peak(
    cluster: list[AdductHypothesis],
    rules: tuple[AdductRule, ...],
) -> dict[int, AdductHypothesis]:
    """Select one deterministic hypothesis per peak for the MSI graph layer.

    CAMERA keeps hypothesis rows until its later checks. The bundled primary
    rules cannot normally assign two rules for one peak within the narrow mass
    tolerance. For an ambiguous custom table, use the hypothesis closest to
    the cluster mean, then CAMERA IPS score, only to obtain one graph label.
    """
    center = float(np.mean([item.neutral_mass for item in cluster]))
    selected: dict[int, AdductHypothesis] = {}
    for item in cluster:
        current = selected.get(item.peak_index)
        if current is None or (
            abs(item.neutral_mass - center),
            -rules[item.rule_index].ips,
            item.rule_index,
        ) < (
            abs(current.neutral_mass - center),
            -rules[current.rule_index].ips,
            current.rule_index,
        ):
            selected[item.peak_index] = item
    return selected


def _candidate_groups_from_clusters(
    clusters: list[list[AdductHypothesis]],
    rules: tuple[AdductRule, ...],
    summed_intensity: np.ndarray,
    *,
    require_quasi: bool,
) -> list[AdductGroup]:
    groups: list[AdductGroup] = []
    for cluster in clusters:
        selected = _choose_rule_per_peak(cluster, rules)
        peak_indices = tuple(sorted(selected))
        rule_indices = tuple(selected[peak].rule_index for peak in peak_indices)
        neutral_masses = tuple(selected[peak].neutral_mass for peak in peak_indices)
        if len(peak_indices) < 2 or len(set(rule_indices)) < 2:
            continue
        if require_quasi and not any(rules[index].quasi for index in rule_indices):
            continue
        representative = _choose_representative_peak(cluster, rules, summed_intensity)
        groups.append(
            AdductGroup(
                neutral_mass=float(np.mean(neutral_masses)),
                peak_indices=peak_indices,
                rule_indices=rule_indices,
                neutral_masses=neutral_masses,
                representative_peak=representative,
                support_peak_indices=peak_indices,
                spatial_correlations=tuple(np.nan for _ in peak_indices),
            )
        )
    return groups


def _validate_shared_mz_data(data_manager: MSDataManagerImzML) -> None:
    if data_manager.ms.meta.continuous is not True:
        raise ValueError(
            "Adduct filtering requires aligned/shared m/z data. "
            "Run peak_align first."
        )
    if data_manager.ms.shared_mz_list is None:
        raise ValueError("Adduct filtering requires a shared m/z axis.")


def _data_manager_state(data_manager: MSDataManagerImzML) -> tuple[object, ...]:
    """Snapshot read-selection state that can change on a manager instance."""
    target_locs = getattr(data_manager, "target_locs", None)
    if target_locs is not None:
        target_locs = tuple(tuple(location) for location in target_locs)
    mz_dtype = getattr(data_manager, "mz_dtype", None)
    intensity_dtype = getattr(data_manager, "intensity_dtype", None)
    return (
        getattr(data_manager, "file_base_path", None),
        target_locs,
        None if mz_dtype is None else np.dtype(mz_dtype).str,
        None if intensity_dtype is None else np.dtype(intensity_dtype).str,
    )


def _candidate_peak_pairs(groups: list[AdductGroup]) -> list[tuple[int, int]]:
    pairs: set[tuple[int, int]] = set()
    for group in groups:
        for left_offset, right_offset in combinations(
            range(len(group.peak_indices)), 2
        ):
            if group.rule_indices[left_offset] == group.rule_indices[right_offset]:
                continue
            left = group.peak_indices[left_offset]
            right = group.peak_indices[right_offset]
            pairs.add((min(left, right), max(left, right)))
    return sorted(pairs)


def _collect_pair_correlations(
    data_manager: MSDataManagerImzML,
    pairs: list[tuple[int, int]],
    *,
    batch_size: int,
) -> dict[tuple[int, int], float]:
    """Stream pairwise Pearson sufficient statistics without storing profiles."""
    if not pairs:
        return {}
    left = np.asarray([pair[0] for pair in pairs], dtype=np.int64)
    right = np.asarray([pair[1] for pair in pairs], dtype=np.int64)
    counts = np.zeros(len(pairs), dtype=np.int64)
    mean_left = np.zeros(len(pairs), dtype=np.float64)
    mean_right = np.zeros(len(pairs), dtype=np.float64)
    m2_left = np.zeros(len(pairs), dtype=np.float64)
    m2_right = np.zeros(len(pairs), dtype=np.float64)
    comoment = np.zeros(len(pairs), dtype=np.float64)

    yielded = False
    for _mz, intensity_flat, lengths, _coords in data_manager.flat_generator(
        batch_size=batch_size,
        include_mz=False,
    ):
        yielded = True
        update_pair_statistics_flat_parallel_jit(
            intensity_flat,
            lengths,
            left,
            right,
            counts,
            mean_left,
            mean_right,
            m2_left,
            m2_right,
            comoment,
        )
    if not yielded:
        raise ValueError("Cannot compute spatial correlations from empty data.")

    denominator = np.sqrt(m2_left * m2_right)
    correlations = np.full(len(pairs), np.nan, dtype=np.float64)
    valid = (counts > 1) & (denominator > 0)
    correlations[valid] = comoment[valid] / denominator[valid]
    np.clip(correlations, -1.0, 1.0, out=correlations)
    return {pair: float(correlations[index]) for index, pair in enumerate(pairs)}


def _maximum_spanning_support(
    component: tuple[int, ...],
    representative: int,
    adjacency: dict[int, dict[int, float]],
) -> dict[int, tuple[int, float]]:
    """Return one accepted correlation edge supporting every component member."""
    support = {representative: (representative, 1.0)}
    visited = {representative}
    remaining = set(component) - visited
    while remaining:
        choices = [
            (correlation, -parent, -child, parent, child)
            for parent in visited
            for child, correlation in adjacency[parent].items()
            if child in remaining
        ]
        if not choices:
            raise RuntimeError(
                "spatial correlation component is unexpectedly disconnected."
            )
        _, _, _, parent, child = max(choices)
        support[child] = (parent, adjacency[parent][child])
        visited.add(child)
        remaining.remove(child)
    return support


def _confirm_groups_from_correlations(
    groups: list[AdductGroup],
    pair_correlations: dict[tuple[int, int], float],
    summed_intensity: np.ndarray,
    rules: tuple[AdductRule, ...],
    *,
    require_quasi: bool,
    min_spatial_correlation: float,
) -> tuple[list[AdductGroup], list[AdductPairEvidence]]:
    """Split mass candidates into MZmine-style correlated components."""
    confirmed: list[AdductGroup] = []
    evidence: list[AdductPairEvidence] = []
    for group in groups:
        rule_by_peak = dict(zip(group.peak_indices, group.rule_indices, strict=True))
        mass_by_peak = dict(zip(group.peak_indices, group.neutral_masses, strict=True))
        adjacency = {peak: {} for peak in group.peak_indices}
        for left, right in combinations(group.peak_indices, 2):
            if rule_by_peak[left] == rule_by_peak[right]:
                continue
            pair = (min(left, right), max(left, right))
            correlation = pair_correlations[pair]
            passes = bool(
                np.isfinite(correlation) and correlation >= min_spatial_correlation
            )
            evidence.append(
                AdductPairEvidence(
                    neutral_mass=group.neutral_mass,
                    left_peak=left,
                    right_peak=right,
                    left_rule=rule_by_peak[left],
                    right_rule=rule_by_peak[right],
                    spatial_correlation=correlation,
                    passes_threshold=passes,
                )
            )
            if passes:
                adjacency[left][right] = correlation
                adjacency[right][left] = correlation

        unvisited = {peak for peak, neighbours in adjacency.items() if neighbours}
        while unvisited:
            start = min(unvisited)
            stack = [start]
            component_set: set[int] = set()
            while stack:
                peak = stack.pop()
                if peak in component_set:
                    continue
                component_set.add(peak)
                stack.extend(adjacency[peak])
            unvisited -= component_set
            component = tuple(sorted(component_set))
            component_rules = tuple(rule_by_peak[peak] for peak in component)
            if len(component) < 2 or len(set(component_rules)) < 2:
                continue
            if require_quasi and not any(
                rules[index].quasi for index in component_rules
            ):
                continue
            representative = min(
                component,
                key=lambda peak: (-float(summed_intensity[peak]), peak),
            )
            support = _maximum_spanning_support(component, representative, adjacency)
            component_masses = tuple(mass_by_peak[peak] for peak in component)
            confirmed.append(
                AdductGroup(
                    neutral_mass=float(np.mean(component_masses)),
                    peak_indices=component,
                    rule_indices=component_rules,
                    neutral_masses=component_masses,
                    representative_peak=representative,
                    support_peak_indices=tuple(support[peak][0] for peak in component),
                    spatial_correlations=tuple(support[peak][1] for peak in component),
                )
            )
    return confirmed, evidence


def _resolve_groups_by_ips(
    groups: list[AdductGroup],
    rules: tuple[AdductRule, ...],
) -> list[AdductGroup]:
    """Apply CAMERA ``checkIps`` competition to overlapping mass groups."""
    unique_groups: list[AdductGroup] = []
    seen: set[tuple[tuple[int, int], ...]] = set()
    for group in groups:
        key = tuple(zip(group.peak_indices, group.rule_indices, strict=True))
        if key not in seen:
            seen.add(key)
            unique_groups.append(group)

    scores = [
        sum(rules[index].ips for index in group.rule_indices) for group in unique_groups
    ]
    peak_to_groups: dict[int, list[int]] = defaultdict(list)
    for group_index, group in enumerate(unique_groups):
        for peak in group.peak_indices:
            peak_to_groups[peak].append(group_index)

    keep = np.ones(len(unique_groups), dtype=bool)
    for group_indices in peak_to_groups.values():
        if len(group_indices) < 2:
            continue
        maximum = max(scores[index] for index in group_indices)
        for index in group_indices:
            if scores[index] < maximum:
                keep[index] = False
    return [group for index, group in enumerate(unique_groups) if keep[index]]


# ---------------------------------------------------------------------------
# Removal mask, report, and public operations
# ---------------------------------------------------------------------------


def _collect_protected_peaks(
    groups: list[AdductGroup] | tuple[AdductGroup, ...],
) -> set[int]:
    return {group.representative_peak for group in groups}


def _build_adduct_keep_mask(
    mz_axis_size: int,
    groups: list[AdductGroup] | tuple[AdductGroup, ...],
) -> np.ndarray:
    protected = _collect_protected_peaks(groups)
    keep_mask = np.ones(mz_axis_size, dtype=bool)
    for group in groups:
        for peak_index in group.peak_indices:
            if peak_index not in protected:
                keep_mask[peak_index] = False
    return keep_mask


def _build_adduct_annotations(
    groups: list[AdductGroup],
    mz_axis: np.ndarray,
    rules: tuple[AdductRule, ...],
    *,
    min_spatial_correlation: float,
) -> list[AdductAnnotation]:
    protected = _collect_protected_peaks(groups)
    annotations: list[AdductAnnotation] = []
    for group in groups:
        rep_offset = group.peak_indices.index(group.representative_peak)
        representative_rule = rules[group.rule_indices[rep_offset]]
        neutral_mass_span = float(max(group.neutral_masses) - min(group.neutral_masses))
        for offset, peak_index in enumerate(group.peak_indices):
            if peak_index == group.representative_peak or peak_index in protected:
                continue
            support_peak = group.support_peak_indices[offset]
            support_offset = group.peak_indices.index(support_peak)
            annotations.append(
                AdductAnnotation(
                    neutral_mass=group.neutral_mass,
                    neutral_mass_span=neutral_mass_span,
                    representative_peak=group.representative_peak,
                    removed_peak=peak_index,
                    representative_mz=float(mz_axis[group.representative_peak]),
                    removed_mz=float(mz_axis[peak_index]),
                    representative_rule=representative_rule.name,
                    removed_rule=rules[group.rule_indices[offset]].name,
                    support_peak=support_peak,
                    support_mz=float(mz_axis[support_peak]),
                    support_rule=rules[group.rule_indices[support_offset]].name,
                    spatial_correlation=float(group.spatial_correlations[offset]),
                    correlation_threshold=min_spatial_correlation,
                )
            )
    return annotations


def _select_summary_annotation_indices(
    annotations: tuple[AdductAnnotation, ...] | list[AdductAnnotation],
) -> set[int]:
    """Select one strongest evidence row for each uniquely removed peak."""
    selected: dict[int, int] = {}
    for index, annotation in enumerate(annotations):
        current_index = selected.get(annotation.removed_peak)
        if current_index is None:
            selected[annotation.removed_peak] = index
            continue
        current = annotations[current_index]
        if (
            -annotation.spatial_correlation,
            annotation.neutral_mass_span,
            annotation.representative_peak,
            annotation.support_peak,
            index,
        ) < (
            -current.spatial_correlation,
            current.neutral_mass_span,
            current.representative_peak,
            current.support_peak,
            current_index,
        ):
            selected[annotation.removed_peak] = index
    return set(selected.values())


def _normalize_adduct_report_path(report_path: str | Path) -> Path:
    path = Path(report_path)
    if not path.suffix:
        return path.with_suffix(".csv")
    if path.suffix.lower() != ".csv":
        raise ValueError("adduct report path must use the .csv extension.")
    return path


def write_adduct_filter_report(
    result: AdductAnnotationResult,
    report_path: str | Path,
) -> Path:
    """Write one four-column summary row per candidate peak."""
    _validate_annotation_result(result)
    path = _normalize_adduct_report_path(report_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "representative_mz",
        "adduct_mz",
        "representative_adduct",
        "adduct",
    ]
    summary_indices = _select_summary_annotation_indices(result.annotations)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for index, annotation in enumerate(result.annotations):
            if index not in summary_indices:
                continue
            writer.writerow(
                {
                    "representative_mz": annotation.representative_mz,
                    "adduct_mz": annotation.removed_mz,
                    "representative_adduct": annotation.representative_rule,
                    "adduct": annotation.removed_rule,
                }
            )
    return path


def _validate_annotation_result(result: AdductAnnotationResult) -> None:
    """Rebuild and validate all removal evidence stored in a preview result."""
    if not isinstance(result, AdductAnnotationResult):
        raise TypeError("result must be an AdductAnnotationResult.")
    if result.mz_axis.ndim != 1:
        raise ValueError("result.mz_axis must be one-dimensional.")
    expected_rule_names = tuple(rule.name for rule in result.rules)
    if result.rule_names != expected_rule_names:
        raise ValueError("result rule names do not match its rule table.")
    n_peaks = result.mz_axis.size
    for group in result.groups:
        lengths = {
            len(group.peak_indices),
            len(group.rule_indices),
            len(group.neutral_masses),
            len(group.support_peak_indices),
            len(group.spatial_correlations),
        }
        if len(lengths) != 1 or group.representative_peak not in group.peak_indices:
            raise ValueError(
                "result contains a structurally inconsistent adduct group."
            )
        if any(not 0 <= peak < n_peaks for peak in group.peak_indices):
            raise ValueError("result contains an out-of-range peak index.")
        if any(not 0 <= rule < len(result.rule_names) for rule in group.rule_indices):
            raise ValueError("result contains an out-of-range rule index.")
        if any(peak not in group.peak_indices for peak in group.support_peak_indices):
            raise ValueError("result contains an invalid support peak.")
        representative_offset = group.peak_indices.index(group.representative_peak)
        if (
            group.support_peak_indices[representative_offset]
            != group.representative_peak
            or group.spatial_correlations[representative_offset] != 1.0
        ):
            raise ValueError("result contains invalid representative evidence.")
        if any(
            not np.isfinite(correlation)
            or correlation < result.config.min_spatial_correlation
            for correlation in group.spatial_correlations
        ):
            raise ValueError("result contains evidence below its recorded threshold.")

    for annotation in result.annotations:
        indices = (
            annotation.representative_peak,
            annotation.removed_peak,
            annotation.support_peak,
        )
        if any(not 0 <= peak < n_peaks for peak in indices):
            raise ValueError("result contains an out-of-range annotation index.")
        if annotation.correlation_threshold != result.config.min_spatial_correlation:
            raise ValueError("result annotation threshold does not match its config.")
        if (
            not np.isfinite(annotation.spatial_correlation)
            or annotation.spatial_correlation < annotation.correlation_threshold
        ):
            raise ValueError("result contains invalid annotation evidence.")

    evidence_keys: set[tuple[float, int, int, int, int]] = set()
    for evidence in result.pair_evidence:
        if not (
            0 <= evidence.left_peak < evidence.right_peak < n_peaks
            and 0 <= evidence.left_rule < len(result.rules)
            and 0 <= evidence.right_rule < len(result.rules)
            and evidence.left_rule != evidence.right_rule
            and np.isfinite(evidence.neutral_mass)
        ):
            raise ValueError("result contains structurally invalid pair evidence.")
        expected_pass = bool(
            np.isfinite(evidence.spatial_correlation)
            and evidence.spatial_correlation >= result.config.min_spatial_correlation
        )
        if evidence.passes_threshold != expected_pass:
            raise ValueError("result pair evidence disagrees with its threshold.")
        evidence_key = (
            evidence.neutral_mass,
            evidence.left_peak,
            evidence.right_peak,
            evidence.left_rule,
            evidence.right_rule,
        )
        if evidence_key in evidence_keys:
            raise ValueError("result contains duplicate pair evidence.")
        evidence_keys.add(evidence_key)

    keep_mask = _build_adduct_keep_mask(n_peaks, result.groups)
    removed_by_mask = set(np.flatnonzero(~keep_mask))
    removed_by_report = {annotation.removed_peak for annotation in result.annotations}
    if removed_by_mask != removed_by_report:
        raise ValueError("result annotations and removal mask are inconsistent.")
    expected_annotations = tuple(
        _build_adduct_annotations(
            list(result.groups),
            result.mz_axis,
            result.rules,
            min_spatial_correlation=result.config.min_spatial_correlation,
        )
    )
    if result.annotations != expected_annotations:
        raise ValueError("result annotation rows do not match its adduct groups.")


def _validate_annotations_belong_to(
    result: AdductAnnotationResult,
    data_manager: MSDataManagerImzML,
) -> None:
    """Validate a preview before applying it to its original data manager."""
    _validate_annotation_result(result)
    if result.source_data_manager is not data_manager:
        raise ValueError(
            "result belongs to a different data manager; recompute the preview "
            "for this exact dataset before applying it."
        )
    if result.source_manager_state != _data_manager_state(data_manager):
        raise ValueError(
            "data manager file, target_locs, or read dtype changed after preview; "
            "recompute the annotations before applying them."
        )
    shared_mz = data_manager.ms.shared_mz_list
    if shared_mz is None or not np.array_equal(result.mz_axis, shared_mz):
        raise ValueError(
            "result.mz_axis does not match the data manager shared m/z axis."
        )


def compute_adduct_annotations(
    data_manager: MSDataManagerImzML,
    *,
    polarity: Polarity = "positive",
    rules: tuple[AdductRule, ...] | None = None,
    ppm: float = 5.0,
    mzabs: float = 0.015,
    require_quasi: bool = True,
    min_spatial_correlation: float = 0.7,
    batch_size: int = 256,
) -> AdductAnnotationResult:
    """Compute auditable adduct decisions without modifying the dataset."""
    _validate_shared_mz_data(data_manager)
    if ppm <= 0 or not np.isfinite(ppm):
        raise ValueError("ppm must be a finite positive number.")
    if mzabs < 0 or not np.isfinite(mzabs):
        raise ValueError("mzabs must be a finite non-negative number.")
    if (
        not np.isfinite(min_spatial_correlation)
        or not -1 <= min_spatial_correlation <= 1
    ):
        raise ValueError("min_spatial_correlation must be finite and in [-1, 1].")
    if not isinstance(batch_size, (int, np.integer)) or batch_size < 1:
        raise ValueError("batch_size must be a positive integer.")

    selected_rules = _resolve_rules(polarity, rules)
    mz_axis = np.asarray(data_manager.ms.shared_mz_list)
    summed_intensity = compute_sum_intensity(data_manager, batch_size=int(batch_size))
    hypotheses = _build_hypotheses(mz_axis, selected_rules)
    clusters = _cluster_hypotheses(hypotheses, ppm=ppm, mzabs=mzabs)
    candidate_groups = _candidate_groups_from_clusters(
        clusters,
        selected_rules,
        summed_intensity,
        require_quasi=require_quasi,
    )
    pairs = _candidate_peak_pairs(candidate_groups)
    pair_correlations = _collect_pair_correlations(
        data_manager,
        pairs,
        batch_size=int(batch_size),
    )
    confirmed_groups, pair_evidence = _confirm_groups_from_correlations(
        candidate_groups,
        pair_correlations,
        summed_intensity,
        selected_rules,
        require_quasi=require_quasi,
        min_spatial_correlation=min_spatial_correlation,
    )
    resolved_groups = _resolve_groups_by_ips(confirmed_groups, selected_rules)
    annotations = _build_adduct_annotations(
        resolved_groups,
        mz_axis,
        selected_rules,
        min_spatial_correlation=min_spatial_correlation,
    )

    immutable_mz_axis = np.array(mz_axis, copy=True)
    immutable_mz_axis.setflags(write=False)
    return AdductAnnotationResult(
        mz_axis=immutable_mz_axis,
        annotations=tuple(annotations),
        groups=tuple(resolved_groups),
        pair_evidence=tuple(pair_evidence),
        rules=selected_rules,
        rule_names=tuple(rule.name for rule in selected_rules),
        config=AdductFilterConfig(
            polarity=polarity,
            ppm=float(ppm),
            mzabs=float(mzabs),
            require_quasi=require_quasi,
            min_spatial_correlation=float(min_spatial_correlation),
            batch_size=int(batch_size),
        ),
        source_data_manager=data_manager,
        source_manager_state=_data_manager_state(data_manager),
    )


def apply_adduct_annotations(
    data_manager: MSDataManagerImzML,
    result: AdductAnnotationResult,
    *,
    output_path: str | Path | None = None,
    temp_dir: str | Path | None = None,
    report_path: str | Path | None = None,
) -> MSDataManagerImzML:
    """Apply a previously previewed result without recomputing evidence."""
    _validate_shared_mz_data(data_manager)
    _validate_annotations_belong_to(result, data_manager)
    if report_path is not None:
        write_adduct_filter_report(result, report_path)
    keep_mask = _build_adduct_keep_mask(result.mz_axis.size, result.groups)
    return filter_data_by_mask(
        data_manager,
        keep_mask,
        output_path=output_path,
        temp_dir=temp_dir,
    )


def filter_peaks_by_adducts(
    data_manager: MSDataManagerImzML,
    *,
    polarity: Polarity = "positive",
    rules: tuple[AdductRule, ...] | None = None,
    ppm: float = 5.0,
    mzabs: float = 0.015,
    require_quasi: bool = True,
    min_spatial_correlation: float = 0.7,
    batch_size: int = 256,
    remove: bool = False,
    output_path: str | Path | None = None,
    temp_dir: str | Path | None = None,
    report_path: str | Path | None = None,
) -> MSDataManagerImzML:
    """Preview/report by default; set ``remove=True`` to apply the decisions."""
    result = compute_adduct_annotations(
        data_manager,
        polarity=polarity,
        rules=rules,
        ppm=ppm,
        mzabs=mzabs,
        require_quasi=require_quasi,
        min_spatial_correlation=min_spatial_correlation,
        batch_size=batch_size,
    )
    if report_path is not None:
        write_adduct_filter_report(result, report_path)
    if not remove:
        return data_manager
    keep_mask = _build_adduct_keep_mask(result.mz_axis.size, result.groups)
    return filter_data_by_mask(
        data_manager,
        keep_mask,
        output_path=output_path,
        temp_dir=temp_dir,
        batch_size=batch_size,
    )
