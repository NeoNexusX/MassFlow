"""Benchmark IonImageBinaryDataset with DataLoader."""

import subprocess
import time
from pathlib import Path

import pytest
from torch.utils.data import DataLoader

from massflow.data_manager import IonImageDataManagerBinary
from massflow.morch import IonImageBinaryDataset

EXAMPLE_IIB = Path("data/example.iib")
ROUND = 1
BATCH_SIZE = 64

_SUDO_PASSWORD = "660789"


def _purge_cache() -> None:
    try:
        subprocess.run(
            ["sudo", "-S", "purge"],
            input=_SUDO_PASSWORD,
            text=True,
            check=True,
            capture_output=True,
            timeout=30,
        )
        time.sleep(1)
    except (subprocess.CalledProcessError, FileNotFoundError):
        pass


@pytest.fixture(scope="class")
def dataset():
    if not EXAMPLE_IIB.exists():
        pytest.skip(f"Binary file not found at {EXAMPLE_IIB}")
    dm = IonImageDataManagerBinary(str(EXAMPLE_IIB), mode="r")
    ds = IonImageBinaryDataset(dm)
    yield ds
    dm.close()


def _run_epoch(dataset, num_workers: int) -> int:
    loader = DataLoader(
        dataset,
        batch_size=BATCH_SIZE,
        shuffle=True,
        num_workers=num_workers,
        pin_memory=False,
    )
    total = 0
    for batch in loader:
        total += batch.size(0)
    return total


class TestIonImageBinaryDatasetBenchmark:
    @pytest.mark.benchmark(timer=time.perf_counter)
    @pytest.mark.parametrize("num_workers", [2])
    def test_dataloader_epoch(self, benchmark, dataset, num_workers):
        total = benchmark.pedantic(
            lambda: _run_epoch(dataset, num_workers),
            rounds=ROUND,
            iterations=1,
            warmup_rounds=1,
            setup=_purge_cache,
        )
        assert total == len(dataset)
