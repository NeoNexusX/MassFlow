import time
import numpy as np
import pytest
from massflow.data_manager import MSDataManagerImzML
from massflow.preprocess.flat_pre_fun import FlatPreprocess
from massflow.tools.logger import get_logger

logger = get_logger("test_normalization")

ROUNDS = 5
FLAT_NORM_METHODS = ["tic_numba", "rms_numba", "ref_numba"]
FLAT_SCALE_NORM_METHODS = ["tic_numba", "rms_numba"]


def _normalization_flat_from_flat_batches(
    flat_batches,
    method: str,
    scale: float | None = None,
    ref_tolerance: float = 0.1,
):
    for mz_flat, intensity_flat, lengths, ref in flat_batches:
        kwargs = {
            "mz_data": mz_flat,
            "intensity": intensity_flat,
            "method": method,
            "scale": scale,
            "lengths": lengths,
        }
        if method == "ref_numba":
            kwargs["mz_flat"] = mz_flat
            kwargs["ref"] = ref
            kwargs["ref_tolerance"] = ref_tolerance

        _ = FlatPreprocess.normalization_flat(**kwargs)


class TestNormalizationAPI:
    """Normalization API tests: memory, speed, consistency, and normalization invariants."""

    @pytest.fixture(scope="module", params=["data/baseline_test.imzML"])
    def flat_caches(self, request):
        data_file_path = request.param
        dm = MSDataManagerImzML(filepath=data_file_path)
        dm.load_head_data()

        caches = []
        for mz_data, intensity_flat, lengths, _ in dm.flat_generator(
            batch_size=4096,
            include_mz=True,
            max_threads=16,
        ):
            ref = float(mz_data[mz_data.size // 2]) if mz_data is not None and mz_data.size > 0 else None
            caches.append((mz_data, intensity_flat, lengths, ref))

        return caches

    @pytest.mark.benchmark(timer=time.perf_counter)
    @pytest.mark.parametrize("method", FLAT_NORM_METHODS)
    def test_norm_flat_speed(
        self,
        benchmark,
        method,
        flat_caches,
    ):
        logger.info("Benchmarking flat normalization method=%s", method)

        benchmark.pedantic(
            _normalization_flat_from_flat_batches,
            args=(flat_caches, method),
            rounds=ROUNDS,
            iterations=1,
            warmup_rounds=1,
        )

    @pytest.mark.parametrize("method", FLAT_SCALE_NORM_METHODS)
    def test_norm_invariants_default_equals_length(self, method, flat_caches):
        _, intensity_flat, lengths, _ = flat_caches[0]
        intensity_flat = np.asarray(intensity_flat, dtype=np.float64)

        normalized_flat = FlatPreprocess.normalization_flat(
            mz_data=None, # type: ignore
            intensity=intensity_flat,
            method=method,
            scale=None,
            lengths=lengths,
        )

        offset = 0
        for valid_len in lengths:
            end = offset + int(valid_len)
            y = normalized_flat.intensity[offset:end]
            if method == "tic_numba":
                assert np.isclose(np.sum(y), float(valid_len), rtol=1e-6)
            else:
                assert np.isclose(np.sqrt(np.mean(y ** 2)), float(valid_len), rtol=1e-6)
            offset = end

    @pytest.mark.parametrize("method", FLAT_SCALE_NORM_METHODS)
    def test_norm_invariants_scale_equals_param(self, method, flat_caches):
        _, intensity_flat, lengths, _ = flat_caches[0]
        intensity_flat = np.asarray(intensity_flat, dtype=np.float64)

        scale_target = 3.5
        normalized_flat = FlatPreprocess.normalization_flat(
            mz_data=None, # type: ignore
            intensity=intensity_flat,
            method=method,
            scale=scale_target,
            lengths=lengths,
        )

        offset = 0
        for valid_len in lengths:
            end = offset + int(valid_len)
            y = normalized_flat.intensity[offset:end]
            if method == "tic_numba":
                assert np.isclose(np.sum(y), scale_target, rtol=1e-6)
            else:
                assert np.isclose(np.sqrt(np.mean(y ** 2)), scale_target, rtol=1e-6)
            offset = end

    def test_norm_ref_peak_scaled_to_target(self, flat_caches):
        mz_flat, intensity_flat, lengths, ref = flat_caches[0]
        assert mz_flat is not None
        assert ref is not None

        is_shared_mz = (
            mz_flat.ndim == 1
            and mz_flat.size != intensity_flat.size
            and lengths.size > 0
            and mz_flat.size == int(lengths[0])
        )

        scale_target = 3.5
        normalized_flat = FlatPreprocess.normalization_flat(
            mz_data=mz_flat,
            intensity=intensity_flat,
            method="ref_numba",
            scale=scale_target,
            mz_flat=mz_flat,
            ref=ref,
            ref_tolerance=0.1,
            lengths=lengths,
        )

        offset = 0
        for valid_len in lengths:
            end = offset + int(valid_len)
            mz_slice = mz_flat if is_shared_mz else mz_flat[offset:end]
            y = normalized_flat.intensity[offset:end]

            idx = int(np.argmin(np.abs(mz_slice - ref)))
            if np.abs(mz_slice[idx] - ref) <= 0.1 and intensity_flat[offset + idx] > 0.0:
                logger.info(
                       f"scaling to target={scale_target},y[idx]={y[idx]}"
                )
                assert np.isclose(y[idx], scale_target, rtol=1e-6, atol=1e-6)

            offset = end
