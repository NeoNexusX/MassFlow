"""Helpers for per-spectrum isotope reporting and optional filtering."""

import csv
from contextlib import ExitStack, contextmanager
from pathlib import Path
from typing import Iterator, Literal, TextIO, cast, overload

import numpy as np

from massflow.data_manager import MSDataManagerImzML
from massflow.preprocess.numba.isotope_filter_numba import (
    apply_isotope_keep_mask_flat_parallel_jit,
    isotope_keep_masks_flat_parallel_jit,
)
from massflow.tools.funs import infer_shared_mz, prepare_flat_inputs

ToleranceUnit = Literal["da", "ppm"]


# ---------------------------------------------------------------------------
# CSV report streaming
# ---------------------------------------------------------------------------


@contextmanager
def _open_isotope_filter_report(
    report_path: str | Path,
) -> Iterator[TextIO]:
    path = Path(report_path)
    if not path.suffix:
        path = path.with_suffix(".csv")
    elif path.suffix.lower() != ".csv":
        raise ValueError("isotope report path must use the .csv extension.")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        csv.writer(handle).writerow(
            ["spectrum_index", "representative_mz", "isotope_mz", "charge"]
        )
        yield handle


def _write_isotope_filter_report_rows(
    handle: TextIO,
    mz_data: np.ndarray,
    lengths: np.ndarray,
    keep: np.ndarray,
    representative_indices: np.ndarray,
    charges: np.ndarray,
    spectrum_offset: int,
) -> None:
    """Append batch rows using spectrum indices relative to the whole input."""
    writer = csv.writer(handle)
    start = 0
    for spectrum_index, length in enumerate(lengths):
        end = start + int(length)
        writer.writerows(
            (
                spectrum_offset + spectrum_index,
                mz_data[start + representative_indices[start + index]],
                mz_data[start + index],
                charges[start + index],
            )
            for index in np.flatnonzero(~keep[start:end])
        )
        start = end


# ---------------------------------------------------------------------------
# Per-spectrum isotope detection and filtering
# ---------------------------------------------------------------------------


def _input_validation(
    mz_arr: np.ndarray,
    lengths_arr: np.ndarray,
    mz_tolerance: float,
    units: ToleranceUnit,
    min_charge: int,
    max_charge: int,
    min_isotope_count: int,
    max_isotope_count: int,
    start_intensity_check: int,
) -> None:
    """Validate inputs for isotope filtering."""
    if infer_shared_mz(mz_arr, lengths_arr):
        raise ValueError("Isotope filtering requires processed + centroided m/z data.")
    if mz_tolerance <= 0 or not np.isfinite(mz_tolerance):
        raise ValueError("mz_tolerance must be a finite positive number.")
    if units not in {"da", "ppm"}:
        raise ValueError("units must be either 'da' or 'ppm'.")
    if (units == "ppm" and mz_tolerance > 100) or (
        units == "da" and mz_tolerance > 0.1
    ):
        raise ValueError("mz_tolerance must not exceed 100 ppm or 0.1 Da.")
    if min_charge < 1 or max_charge < min_charge:
        raise ValueError("Charges must satisfy 1 <= min_charge <= max_charge.")
    if (
        min_isotope_count < 2
        or max_isotope_count < min_isotope_count
    ):
        raise ValueError(
            "Isotope counts must satisfy 2 <= min_isotope_count <= max_isotope_count."
        )
    if start_intensity_check < 0:
        raise ValueError("start_intensity_check must be non-negative.")


def _filter_isotope_batch(
    mz_data: np.ndarray,
    intensity: np.ndarray,
    lengths: np.ndarray,
    *,
    mz_tolerance: float = 1.0,
    units: Literal["da", "ppm"] = "ppm",
    min_charge: int = 1,
    max_charge: int = 3,
    min_isotope_count: int = 3,
    max_isotope_count: int = 10,
    use_decreasing_model: bool = True,
    start_intensity_check: int = 2,
    remove: bool = False,
    report_handle: TextIO | None = None,
    spectrum_offset: int = 0,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Compute one batch once, then report and optionally compact its arrays."""
    mz_arr, intensity_arr, lengths_arr = prepare_flat_inputs(
        mz_data,
        intensity,
        lengths,
    )
    _input_validation(
        mz_arr,
        lengths_arr,
        mz_tolerance,
        units,
        min_charge,
        max_charge,
        min_isotope_count,
        max_isotope_count,
        start_intensity_check
    )

    keep, filtered_lengths, representative_indices, charges = (
        isotope_keep_masks_flat_parallel_jit(
            mz_arr,
            intensity_arr,
            lengths_arr,
            mz_tolerance,
            units == "ppm",
            min_charge,
            max_charge,
            min_isotope_count,
            max_isotope_count,
            use_decreasing_model,
            start_intensity_check,
            report_handle is not None,
        )
    )
    if report_handle is not None:
        _write_isotope_filter_report_rows(
            report_handle,
            mz_arr,
            lengths_arr,
            keep,
            representative_indices,
            charges,
            spectrum_offset,
        )
    if not remove:
        return mz_arr, intensity_arr, lengths_arr
    filtered_mz, filtered_intensity = apply_isotope_keep_mask_flat_parallel_jit(
        mz_arr,
        intensity_arr,
        lengths_arr,
        keep,
        filtered_lengths,
    )
    return filtered_mz, filtered_intensity, filtered_lengths


# ---------------------------------------------------------------------------
# Standalone array and imzML entry point
# ---------------------------------------------------------------------------


@overload
def filter_peaks_by_isotopes(
    mz_data: np.ndarray,
    intensity: np.ndarray,
    lengths: np.ndarray,
    *,
    mz_tolerance: float = 1.0,
    units: ToleranceUnit = "ppm",
    min_charge: int = 1,
    max_charge: int = 3,
    min_isotope_count: int = 3,
    max_isotope_count: int = 10,
    use_decreasing_model: bool = True,
    start_intensity_check: int = 2,
    remove: bool = False,
    report_path: str | Path | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]: ...


@overload
def filter_peaks_by_isotopes(
    *,
    data_manager: MSDataManagerImzML,
    mz_tolerance: float = 1.0,
    units: ToleranceUnit = "ppm",
    min_charge: int = 1,
    max_charge: int = 3,
    min_isotope_count: int = 3,
    max_isotope_count: int = 10,
    use_decreasing_model: bool = True,
    start_intensity_check: int = 2,
    remove: bool = False,
    report_path: str | Path | None = None,
    output_path: str | Path | None = None,
    temp_dir: str | Path | None = None,
    batch_size: int = 256,
) -> MSDataManagerImzML: ...


def filter_peaks_by_isotopes(
    mz_data: np.ndarray | None = None,
    intensity: np.ndarray | None = None,
    lengths: np.ndarray | None = None,
    *,
    data_manager: MSDataManagerImzML | None = None,
    mz_tolerance: float = 1.0,
    units: ToleranceUnit = "ppm",
    min_charge: int = 1,
    max_charge: int = 3,
    min_isotope_count: int = 3,
    max_isotope_count: int = 10,
    use_decreasing_model: bool = True,
    start_intensity_check: int = 2,
    remove: bool = False,
    report_path: str | Path | None = None,
    output_path: str | Path | None = None,
    temp_dir: str | Path | None = None,
    batch_size: int = 256,
) -> MSDataManagerImzML | tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Report isotope peaks from ``data_manager=dm`` or three flat input arrays.

    The imzML input must have loaded processed + centroided metadata, before
    peak alignment. Its spectra and report rows are streamed by ``batch_size``.
    By default the original manager is returned without creating an output
    dataset. ``remove=True`` returns a new manager written to ``output_path``
    (or temporary storage when omitted). The caller owns the returned manager.

    Array input returns the three arrays, filtered only with ``remove=True``.
    Clusters use successive C13-C12 mass differences divided by charge, keeping
    the lowest-m/z peak and all unassigned peaks when filtering. ``report_path``
    writes a four-column CSV with zero-based spectrum indices across batches.
    """
    if data_manager is not None:
        meta = data_manager.ms.meta
        if meta.continuous is True or meta.centroid_spectrum is not True:
            raise ValueError(
                "Isotope filtering requires processed + centroided imzML data; "
                "load metadata and perform peak picking before this step."
            )

    with ExitStack() as stack:
        report_handle = (
            stack.enter_context(_open_isotope_filter_report(report_path))
            if report_path is not None else None
        )

        def filter_batch(
            batch_mz: np.ndarray,
            batch_intensity: np.ndarray,
            batch_lengths: np.ndarray,
            spectrum_offset: int = 0,
        ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
            return _filter_isotope_batch(
                batch_mz, batch_intensity, batch_lengths,
                mz_tolerance=mz_tolerance,
                units=units,
                min_charge=min_charge,
                max_charge=max_charge,
                min_isotope_count=min_isotope_count,
                max_isotope_count=max_isotope_count,
                use_decreasing_model=use_decreasing_model,
                start_intensity_check=start_intensity_check,
                remove=remove,
                report_handle=report_handle,
                spectrum_offset=spectrum_offset,
            )

        if data_manager is None:
            return filter_batch(
                cast(np.ndarray, mz_data),
                cast(np.ndarray, intensity),
                cast(np.ndarray, lengths),
            )

        result = data_manager
        if remove:
            if output_path is not None:
                Path(output_path).parent.mkdir(parents=True, exist_ok=True)
            result = MSDataManagerImzML(
                filepath=None if output_path is None else str(output_path),
                temp_dir=None if temp_dir is None else str(temp_dir),
                mz_dtype=data_manager.mz_dtype,
                intensity_dtype=data_manager.intensity_dtype,
            )
            if (
                Path(result.imzml_filepath).resolve()
                == Path(data_manager.imzml_filepath).resolve()
            ):
                raise ValueError("output_path must differ from the input imzML path.")
            result.copy_meta(data_manager)
            result.ms.meta.processed = True
            result.ms.meta.continuous = None
            result.ms.meta.centroid_spectrum = True
            result.ms.meta.profile_spectrum = None
            stack.callback(result.close_writer)

        spectrum_offset = 0
        for (
            batch_mz, batch_intensity, batch_lengths, coordinates,
        ) in data_manager.flat_generator(
            batch_size=batch_size, include_mz=True,
        ):
            filtered_mz, filtered_intensity, filtered_lengths = filter_batch(
                batch_mz, batch_intensity, batch_lengths, spectrum_offset,
            )
            if remove:
                result.write_flat_batch(
                    intensity_flat=filtered_intensity,
                    lengths=filtered_lengths,
                    coordinates=coordinates,
                    mz_axis=filtered_mz,
                )
            spectrum_offset += batch_lengths.size
    if remove:
        result.load_head_data()
    return result
