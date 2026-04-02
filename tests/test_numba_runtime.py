import pytest
from numba import get_num_threads, set_num_threads

from massflow.preprocess.numba.numba_runtime import (
    apply_numba_runtime,
    detect_performance_core_workers,
    get_global_numba_runtime,
    get_logical_cpu_count,
    set_global_numba_runtime,
)


@pytest.fixture(autouse=True)
def _restore_numba_runtime_state():
    prev_runtime_threads = get_num_threads()
    prev_global_runtime = get_global_numba_runtime()
    yield
    set_global_numba_runtime(prev_global_runtime["max_workers"])
    set_num_threads(prev_runtime_threads)


def test_set_global_numba_runtime_roundtrip() -> None:
    target = 1

    configured = set_global_numba_runtime(target)

    assert configured["max_workers"] == target
    assert get_global_numba_runtime()["max_workers"] == target
    assert get_num_threads() == target


def test_apply_numba_runtime_uses_global_when_not_provided() -> None:
    target = 1

    set_global_numba_runtime(target)
    applied = apply_numba_runtime()

    assert applied == target
    assert get_num_threads() == target


def test_apply_numba_runtime_prefers_explicit_override() -> None:
    prev_threads = get_num_threads()
    supports_two_threads = True
    try:
        set_num_threads(2)
    except ValueError:
        supports_two_threads = False
    finally:
        set_num_threads(prev_threads)

    if not supports_two_threads:
        pytest.skip("Need at least 2 threads to verify explicit override precedence")

    set_global_numba_runtime(1)
    applied = apply_numba_runtime(2)

    assert applied == 2
    assert get_num_threads() == 2


def test_set_global_numba_runtime_validates_range() -> None:
    with pytest.raises(ValueError):
        set_global_numba_runtime(0)

    with pytest.raises(ValueError):
        set_global_numba_runtime(10**9)


def test_set_global_numba_runtime_validates_type() -> None:
    with pytest.raises(TypeError):
        set_global_numba_runtime(1.5)  # type: ignore[arg-type]


def test_get_logical_cpu_count_returns_positive_int() -> None:
    count = get_logical_cpu_count()
    print(f"Logical CPU count: {count}")
    assert isinstance(count, int)
    assert count >= 1


def test_detect_performance_core_workers_returns_positive_int() -> None:
    count = detect_performance_core_workers()
    assert isinstance(count, int)
    assert count >= 1
