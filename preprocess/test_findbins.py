"""
Test script for `_findbins` function in `preprocess/est_noise.py`.

This script provides multiple scenarios to verify that `_findbins` behaves
correctly under different configurations (static binning with and without
overlap, dynamic binning, limits-only mode, invalid input handling, and
nbins clipping). It prints PASS/FAIL per test and sets a non-zero exit code
on failure.

Run: `python preprocess/test_findbins.py`
"""

from typing import Callable, Tuple
import sys
import numpy as np

from preprocess.est_noise import _findbins


def generate_signal(n: int) -> np.ndarray:
    """Generate a simple 1D synthetic signal for testing.

    The signal combines a linear trend and a small sinusoidal variation, which
    provides non-trivial structure for binning and SSE calculations.

    Args:
        n (int): Length of the signal. Must be a positive integer.

    Returns:
        np.ndarray: 1D array of length `n` representing the synthetic signal.

    Raises:
        ValueError: If `n` is not a positive integer.
    """
    if not isinstance(n, int) or n <= 0:
        raise ValueError("n must be a positive integer")
    x = np.arange(n, dtype=float)
    return x + 0.1 * np.sin(2 * np.pi * x / max(n, 1))


def run_test(name: str, func: Callable[[], None]) -> bool:
    """Run a single test function and print PASS/FAIL.

    Args:
        name (str): Human-readable test name.
        func (Callable[[], None]): Test callable that raises on failure.

    Returns:
        bool: True if the test passes without exception, False otherwise.

    Raises:
        None: This function catches all exceptions from `func` and does not propagate.
    """
    try:
        func()
        print(f"[PASS] {name}")
        return True
    except Exception as e:
        print(f"[FAIL] {name}: {e}")
        return False


def test_static_no_overlap() -> None:
    """Validate static binning with no overlap produces contiguous partitions.

    Scenario:
    - Use `overlap=0` and a small signal; expect `nbins` disjoint bins that
      cover the entire index range without gaps or overlaps.

    Assertions:
    - Number of bins equals `nbins`.
    - Lower/upper boundaries cover [0, n-1] exactly.
    - Contiguity: `lower[i+1] == upper[i] + 1`.

    Args:
        None

    Returns:
        None

    Raises:
        AssertionError: If any of the assertions fail.
    """
    n = 10
    nbins = 3
    data = generate_signal(n)
    bins, meta = _findbins(data, nbins=nbins, overlap=0, dynamic=False)

    lower = meta["lower"]
    upper = meta["upper"]
    sizes = meta["size"]

    assert len(bins) == nbins, "Unexpected number of bins"
    assert lower[0] == 0 and upper[-1] == n - 1, "Coverage mismatch"
    assert np.all(lower[1:] == upper[:-1] + 1), "Bins are not contiguous"
    # Bin slices match sizes and data slices
    for i, (li, ui) in enumerate(zip(lower, upper)):
        expected = data[int(li): int(ui) + 1]
        assert len(bins[i]) == sizes[i], "Size metadata mismatch"
        assert np.allclose(bins[i], expected), "Bin slice content mismatch"


def test_static_with_overlap() -> None:
    """Validate static binning with overlap produces overlapping bins correctly.

    Scenario:
    - Use `overlap=0.5` with `dynamic=False`; expect `nbins` bins whose sizes
      may sum to more than `n` and include overlaps.

    Assertions:
    - Number of bins equals `nbins`.
    - Every bin has `upper[i] >= lower[i]` and size >= 1.
    - Bin slices match `meta` boundaries.

    Args:
        None

    Returns:
        None

    Raises:
        AssertionError: If any of the assertions fail.
    """
    n = 10
    nbins = 2
    data = generate_signal(n)
    bins, meta = _findbins(data, nbins=nbins, overlap=0.5, dynamic=False)

    lower = meta["lower"]
    upper = meta["upper"]
    sizes = meta["size"]

    assert len(bins) == nbins, "Unexpected number of bins"
    for i, (li, ui) in enumerate(zip(lower, upper)):
        assert ui >= li, "Upper bound should be >= lower bound"
        assert sizes[i] >= 1, "Bin size must be positive"
        expected = data[int(li): int(ui) + 1]
        assert len(bins[i]) == sizes[i], "Size metadata mismatch"
        assert np.allclose(bins[i], expected), "Bin slice content mismatch"


def test_dynamic_partitioning() -> None:
    """Validate dynamic binning partitions the signal and returns optimization info.

    Scenario:
    - Use `dynamic=True` with `nbins` intentionally < 3 to trigger clip to 3.
    - Expect contiguous coverage across bins and presence of `sse` and `trace`.

    Assertions:
    - Number of bins equals clipped `nbins` (>=3).
    - Coverage `[0, n-1]` and contiguity across bins.
    - `sse` and `trace` exist; `trace` is non-increasing overall.

    Args:
        None

    Returns:
        None

    Raises:
        AssertionError: If any of the assertions fail.
    """
    n = 30
    nbins_requested = 2  # will be clipped to 3 in dynamic mode
    data = generate_signal(n)
    bins, meta = _findbins(data, nbins=nbins_requested, dynamic=True, niter=10)

    lower = meta["lower"]
    upper = meta["upper"]
    sizes = meta["size"]
    sse = meta.get("sse")
    trace = meta.get("trace")

    # nbins clipped to 3
    assert len(bins) == 3, "Dynamic mode should clip nbins to >= 3"
    assert lower[0] == 0 and upper[-1] == n - 1, "Coverage mismatch"
    assert np.all(lower[1:] == upper[:-1] + 1), "Bins are not contiguous"
    assert sse is not None and trace is not None, "Missing dynamic extras"
    assert len(trace) >= 1, "Optimization trace should have at least one entry"
    # Non-increasing trace (allow equal due to acceptance criteria)
    assert trace[0] >= trace[-1], "Trace did not improve or stay the same"
    # Validate slices match meta boundaries
    for i, (li, ui) in enumerate(zip(lower, upper)):
        expected = data[int(li): int(ui) + 1]
        assert len(bins[i]) == sizes[i], "Size metadata mismatch"
        assert np.allclose(bins[i], expected), "Bin slice content mismatch"


def test_limits_only() -> None:
    """Validate `limits_only=True` returns only metadata.

    Scenario:
    - Request limits-only output; expect a dict containing lower/upper/size and
      no bins list.

    Assertions:
    - Return type is `dict`.
    - Contains keys: `lower`, `upper`, `size`.

    Args:
        None

    Returns:
        None

    Raises:
        AssertionError: If any of the assertions fail.
    """
    n = 12
    data = generate_signal(n)
    meta = _findbins(data, nbins=4, overlap=0, dynamic=False, limits_only=True)
    assert isinstance(meta, dict), "limits_only should return a metadata dict"
    for key in ("lower", "upper", "size"):
        assert key in meta, f"Missing key in metadata: {key}"


def test_invalid_ndim_raises() -> None:
    """Validate invalid data dimensionality raises `ValueError`.

    Scenario:
    - Pass a 2D array; expect `ValueError` indicating data must be 1D.

    Args:
        None

    Returns:
        None

    Raises:
        AssertionError: If a `ValueError` is not raised.
    """
    bad = np.zeros((5, 5))
    try:
        _findbins(bad, nbins=2, overlap=0)
    except ValueError:
        return
    raise AssertionError("Expected ValueError for non-1D input")


def test_nbins_clip_to_length() -> None:
    """Validate `nbins` is clipped to data length when larger than `n`.

    Scenario:
    - Use `nbins` much larger than `n` with `overlap=0`; expect `nbins == n`.

    Args:
        None

    Returns:
        None

    Raises:
        AssertionError: If `nbins` is not clipped appropriately.
    """
    n = 5
    nbins_requested = 100
    data = generate_signal(n)
    bins, meta = _findbins(data, nbins=nbins_requested, overlap=0, dynamic=False)
    assert len(bins) == n, "nbins should be clipped to data length"
    lower = meta["lower"]
    upper = meta["upper"]
    # Each bin should be a single element partition for n bins over n points
    assert np.all(upper - lower + 1 == 1), "Expected unit-size bins after clipping"


def main() -> None:
    """Run all `_findbins` tests and exit with appropriate status code.

    Args:
        None

    Returns:
        None

    Raises:
        SystemExit: Exits with status code 0 if all tests pass, 1 otherwise.
    """
    tests: Tuple[Tuple[str, Callable[[], None]], ...] = (
        ("static_no_overlap", test_static_no_overlap),
        ("static_with_overlap", test_static_with_overlap),
        ("dynamic_partitioning", test_dynamic_partitioning),
        ("limits_only", test_limits_only),
        ("invalid_ndim_raises", test_invalid_ndim_raises),
        ("nbins_clip_to_length", test_nbins_clip_to_length),
    )

    results = [run_test(name, func) for name, func in tests]
    all_ok = all(results)
    if all_ok:
        print("All tests passed.")
        raise SystemExit(0)
    else:
        print("Some tests failed.")
        raise SystemExit(1)


if __name__ == "__main__":
    main()