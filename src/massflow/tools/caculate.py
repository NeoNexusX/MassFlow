import os
from pathlib import Path
import argparse
import numpy as np
from pyimzml.ImzMLParser import ImzMLParser

DEFAULT_IMZML = ""
DEFAULT_SCAN_ALL = True
DEFAULT_SAMPLE = 20

def _dtype_bytes(val):
    if val is None:
        return None
    if isinstance(val, np.dtype):
        return val.itemsize
    if isinstance(val, type) and issubclass(val, np.generic):
        return np.dtype(val).itemsize
    if isinstance(val, int):
        if val in (8, 16, 32, 64):
            return val // 8
        return None
    if isinstance(val, str):
        try:
            return np.dtype(val).itemsize
        except Exception:
            return None
    return None

def _format_size(nbytes):
    if nbytes is None:
        return "unknown"
    units = ["B", "KiB", "MiB", "GiB", "TiB"]
    size = float(nbytes)
    idx = 0
    while size >= 1024.0 and idx < len(units) - 1:
        size /= 1024.0
        idx += 1
    return f"{size:.3f} {units[idx]}"

def _get_lengths(parser, scan_all, sample, spectra_count):
    if spectra_count == 0:
        return 0, 0, False
    mz0, _ = parser.getspectrum(0)
    len0 = len(mz0)
    if scan_all:
        total_points = 0
        for i in range(spectra_count):
            mz, _ = parser.getspectrum(i)
            total_points += len(mz)
        return len0, total_points, True
    sample_n = min(sample, spectra_count)
    if sample_n <= 1:
        return len0, len0 * spectra_count, True
    step = max(1, spectra_count // sample_n)
    same = True
    for i in range(0, spectra_count, step):
        mz, _ = parser.getspectrum(i)
        if len(mz) != len0:
            same = False
            break
    if not same:
        return len0, len0 * spectra_count, False
    return len0, len0 * spectra_count, True

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--imzml", default=DEFAULT_IMZML)
    ap.add_argument("--mode", choices=["processed", "continuous"], default="processed")
    ap.add_argument("--scan-all", action="store_true", default=DEFAULT_SCAN_ALL)
    ap.add_argument("--sample", type=int, default=DEFAULT_SAMPLE)
    args = ap.parse_args()

    imzml_path = Path(args.imzml)
    if not imzml_path.exists():
        raise FileNotFoundError(f"imzML not found: {imzml_path}")

    parser = ImzMLParser(str(imzml_path))
    spectra_count = len(parser.coordinates)

    mz_bytes = _dtype_bytes(getattr(parser, "mzPrecision", None))
    inten_bytes = _dtype_bytes(getattr(parser, "intensityPrecision", None))

    if mz_bytes is None or inten_bytes is None:
        mz0, inten0 = parser.getspectrum(0)
        mz_bytes = mz_bytes or mz0.dtype.itemsize
        inten_bytes = inten_bytes or inten0.dtype.itemsize

    len0, total_points, length_consistent = _get_lengths(
        parser, args.scan_all, args.sample, spectra_count
    )

    if args.mode == "processed":
        total_mz_bytes = total_points * mz_bytes
    else:
        total_mz_bytes = len0 * mz_bytes

    total_int_bytes = total_points * inten_bytes
    total_bytes = total_mz_bytes + total_int_bytes

    imzml_size = imzml_path.stat().st_size
    ibd_path = imzml_path.with_suffix(".ibd")
    ibd_size = ibd_path.stat().st_size if ibd_path.exists() else None

    print(f"imzML: {imzml_path}")
    print(f"spectra: {spectra_count}")
    print(f"mode: {args.mode}")
    print(f"mz dtype bytes: {mz_bytes}")
    print(f"intensity dtype bytes: {inten_bytes}")
    print(f"spectrum length (first): {len0}")
    print(f"total points: {total_points}")
    print(f"length consistent: {length_consistent}")
    if not length_consistent and not args.scan_all:
        print("warning: variable spectrum lengths detected in samples, use --scan-all for exact total points")
    print(f"raw mz bytes: {total_mz_bytes} ({_format_size(total_mz_bytes)})")
    print(f"raw intensity bytes: {total_int_bytes} ({_format_size(total_int_bytes)})")
    print(f"raw total bytes: {total_bytes} ({_format_size(total_bytes)})")
    print(f"imzML file size: {imzml_size} ({_format_size(imzml_size)})")
    if ibd_size is not None:
        print(f"ibd file size: {ibd_size} ({_format_size(ibd_size)})")
        print(f"imzML+ibd size: {imzml_size + ibd_size} ({_format_size(imzml_size + ibd_size)})")

if __name__ == "__main__":
    main()