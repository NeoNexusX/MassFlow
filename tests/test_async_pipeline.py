import time
from pathlib import Path
import numpy as np
import pytest
from massflow.data_manager import (
    MSDataManagerImzML,
    IonImageDataManagerBinary,
    IonImageDataManagerH5,
    IonImageDataManagerZarr,
)
from massflow.preprocess import Preprocessor
from massflow.tools.logger import get_logger

logger = get_logger("massflow.tests.test_async_pipeline")
ROUND = 1
ION_NUM = 32


def _run_pipeline(
    data_manager,
    batch_size=256,
    queue_ab_size=2,
    queue_bc_size=2,
    output_format="imzml",
    temp_dir="./temp",
    output_path=None,
    compute_mean_spectrum=False,
):
    """Run the async pipeline once with specified parameters."""
    processed_data_manager = (
        Preprocessor(
            data_manager=data_manager,
            batch_size=batch_size,
            queue_ab_size=queue_ab_size,
            queue_bc_size=queue_bc_size,
            temp_dir=str(temp_dir),
            output_format=output_format, # type: ignore
            output_path=None if output_path is None else str(output_path),
            compute_mean_spectrum=compute_mean_spectrum,
        )
        .noise_reduction()
        .peak_pick()
        .peak_align()
        .start()
    )
    return processed_data_manager

def _del_processed_dm_file(processed_dm):
    """Helper to delete the output file of the processed data manager if it exists."""
    if isinstance(processed_dm, MSDataManagerImzML):
        processed_dm.close()
    else:
        processed_dm.finalize()

def _close_after_test(request, data_manager):
    request.addfinalizer(data_manager.close)
    return data_manager

def _sample_indices(size: int, sample_count: int = ION_NUM) -> list[int]:
    sample_count = min(size, sample_count)
    if sample_count <= 0:
        return []
    if sample_count == 1:
        return [0]
    return np.linspace(0, size - 1, sample_count, dtype=np.int64).tolist()

def _collect_imzml_matrix(data_manager: MSDataManagerImzML):
    mz_axis = None
    pixel_batches = []

    for batch_mz, intensity_flat, lengths, _coords in data_manager.flat_generator(
        batch_size=128,
        include_mz=True,
        max_threads=2,
    ):
        batch_mz = np.asarray(batch_mz)
        mz_axis = batch_mz.copy() if mz_axis is None else mz_axis
        np.testing.assert_allclose(batch_mz, mz_axis, rtol=1e-7, atol=1e-10)

        n_ions = len(mz_axis)
        assert np.all(np.asarray(lengths) == n_ions)
        pixel_batches.append(np.asarray(intensity_flat).reshape(len(lengths), n_ions))

    assert mz_axis is not None
    return mz_axis, np.vstack(pixel_batches)

def _assert_ion_image_matches(ion_dm, mz_axis, pixel_matrix):
    assert ion_dm.n_ions == len(mz_axis)
    assert ion_dm.n_pixels == pixel_matrix.shape[0]
    assert ion_dm.n_ions > 0
    assert ion_dm.n_pixels > 0
    np.testing.assert_allclose(ion_dm.mz_axis, mz_axis, rtol=1e-7, atol=1e-10)

    for ion_idx in _sample_indices(ion_dm.n_ions):
        np.testing.assert_allclose(
            ion_dm.get_ion_vector(ion_idx),
            pixel_matrix[:, ion_idx],
            rtol=1e-5,
            atol=1e-8,
        )

class TestAsyncPipeline:
    """Test class for the massflow pipeline."""

    @pytest.fixture(scope="module", params=["data/example.imzML"])
    def ms_raw_data(self, request) -> MSDataManagerImzML:
        """Fixture providing MSDataManagerImzML instance with fully initialized spectra for noise reduction tests."""
        data_file_path = Path(request.param)
        if not data_file_path.exists() or not data_file_path.with_suffix(".ibd").exists():
            pytest.skip(f"Missing test data: {data_file_path} and matching .ibd file")

        dm = MSDataManagerImzML(filepath=data_file_path)
        dm.load_head_data()
        request.addfinalizer(dm.close)
        return dm

    @pytest.mark.parametrize("batch_size", [256, 512])
    @pytest.mark.parametrize("queue_ab_size", [2])
    @pytest.mark.parametrize("queue_bc_size", [2])
    @pytest.mark.parametrize("output_format", ["imzml", "ion_image_binary", "ion_image_h5", "ion_image_zarr"])
    @pytest.mark.benchmark(timer=time.perf_counter)
    def test_async_pipeline_benchmark(
        self,
        ms_raw_data,
        benchmark,
        batch_size: int,
        queue_ab_size: int,
        queue_bc_size: int,
        output_format: str,
    ) -> None:
        """Whole pipeline runtime benchmark."""
        processed_dm = benchmark.pedantic(
            _run_pipeline,
            args=(ms_raw_data, batch_size, queue_ab_size, queue_bc_size, output_format),
            rounds=ROUND,
            iterations=1,
            warmup_rounds=1,
        )
        _del_processed_dm_file(processed_dm)

    def test_binary_matches_imzml(self, ms_raw_data, tmp_path, request):
        """Binary output must match the same pipeline's imzML output."""
        imzml_dm = _close_after_test(
            request,
            _run_pipeline(ms_raw_data, temp_dir=tmp_path / "imzml"),
        )
        mz_axis, pixel_matrix = _collect_imzml_matrix(imzml_dm) # type: ignore

        output_path = tmp_path / "pipeline.iib"
        binary_dm = _close_after_test(
            request,
            _run_pipeline(
                ms_raw_data,
                output_format="ion_image_binary",
                temp_dir=tmp_path / "binary",
                output_path=output_path,
            ),
        )

        assert isinstance(binary_dm, IonImageDataManagerBinary)
        assert output_path.exists()
        _assert_ion_image_matches(binary_dm, mz_axis, pixel_matrix)
        for pixel_idx in _sample_indices(binary_dm.n_pixels):
            np.testing.assert_allclose(
                binary_dm.get_pixel_spectrum(pixel_idx),
                pixel_matrix[pixel_idx],
                rtol=1e-5,
                atol=1e-8,
            )

    def test_h5_matches_imzml(self, ms_raw_data, tmp_path, request):
        """H5 output must match the same pipeline's imzML output."""
        imzml_dm = _close_after_test(
            request,
            _run_pipeline(ms_raw_data, temp_dir=tmp_path / "imzml"),
        )
        mz_axis, pixel_matrix = _collect_imzml_matrix(imzml_dm) # type: ignore

        output_path = tmp_path / "pipeline.h5"
        h5_dm = _close_after_test(
            request,
            _run_pipeline(
                ms_raw_data,
                output_format="ion_image_h5",
                temp_dir=tmp_path / "h5",
                output_path=output_path,
            ),
        )

        assert isinstance(h5_dm, IonImageDataManagerH5)
        assert output_path.exists()
        _assert_ion_image_matches(h5_dm, mz_axis, pixel_matrix)

    def test_zarr_matches_imzml(self, ms_raw_data, tmp_path, request):
        """zarr output must match the same pipeline's imzML output."""
        imzml_dm = _close_after_test(
            request,
            _run_pipeline(ms_raw_data, temp_dir=tmp_path / "imzml"),
        )
        mz_axis, pixel_matrix = _collect_imzml_matrix(imzml_dm) # type: ignore

        output_path = tmp_path / "pipeline.zarr"
        zarr_dm = _close_after_test(
            request,
            _run_pipeline(
                ms_raw_data,
                output_format="ion_image_zarr",
                temp_dir=tmp_path / "zarr",
                output_path=output_path,
            ),
        )

        assert isinstance(zarr_dm, IonImageDataManagerZarr)
        assert output_path.exists()
        _assert_ion_image_matches(zarr_dm, mz_axis, pixel_matrix)

    def test_zarr_pipeline_computes_mean_spectrum(self, ms_raw_data, tmp_path, request):
        output_path = tmp_path / "pipeline_mean.zarr"
        zarr_dm = _close_after_test(
            request,
            _run_pipeline(
                ms_raw_data,
                output_format="ion_image_zarr",
                temp_dir=tmp_path / "zarr_mean",
                output_path=output_path,
                compute_mean_spectrum=True,
            ),
        )

        assert isinstance(zarr_dm, IonImageDataManagerZarr)
        ion_indices = _sample_indices(zarr_dm.n_ions)
        expected = np.array(
            [zarr_dm.get_ion_vector(ion_idx).mean(dtype=np.float64) for ion_idx in ion_indices]
        )
        np.testing.assert_allclose(zarr_dm.mean_spectrum[ion_indices], expected, rtol=1e-7, atol=1e-10) # type: ignore

    def test_zarr_conversion_computes_mean_spectrum(self, ms_raw_data, tmp_path, request):
        output_path = tmp_path / "conversion_mean.zarr"
        zarr_dm = _close_after_test(
            request,
            Preprocessor(
                data_manager=ms_raw_data,
                temp_dir=str(tmp_path / "conversion_mean"),
                output_format="ion_image_zarr",
                output_path=str(output_path),
                compute_mean_spectrum=True,
            ).start(),
        )

        assert isinstance(zarr_dm, IonImageDataManagerZarr)
        ion_indices = _sample_indices(zarr_dm.n_ions)
        expected = np.array(
            [zarr_dm.get_ion_vector(ion_idx).mean(dtype=np.float64) for ion_idx in ion_indices]
        )
        np.testing.assert_allclose(zarr_dm.mean_spectrum[ion_indices], expected, rtol=1e-7, atol=1e-10) # type: ignore

    def test_three_output_formats_are_consistent(self, ms_raw_data, tmp_path, request):
        """The three output formats must be consistent with each other."""
        binary_path = tmp_path / "pipeline.iib"
        binary_dm = _close_after_test(
            request,
            _run_pipeline(
                ms_raw_data,
                output_format="ion_image_binary",
                temp_dir=tmp_path / "binary",
                output_path=binary_path,
            ),
        )

        h5_path = tmp_path / "pipeline.h5"
        h5_dm = _close_after_test(
            request,
            _run_pipeline(
                ms_raw_data,
                output_format="ion_image_h5",
                temp_dir=tmp_path / "h5",
                output_path=h5_path,
            ),
        )

        zarr_path = tmp_path / "pipeline.zarr"
        zarr_dm = _close_after_test(
            request,
            _run_pipeline(
                ms_raw_data,
                output_format="ion_image_zarr",
                temp_dir=tmp_path / "zarr",
                output_path=zarr_path,
            ),
        )

        assert isinstance(binary_dm, IonImageDataManagerBinary)
        assert isinstance(h5_dm, IonImageDataManagerH5)
        assert isinstance(zarr_dm, IonImageDataManagerZarr)
        assert binary_path.exists()
        assert h5_path.exists()
        assert zarr_path.exists()

        assert h5_dm.n_ions == binary_dm.n_ions == zarr_dm.n_ions
        assert h5_dm.n_pixels == binary_dm.n_pixels == zarr_dm.n_pixels
        np.testing.assert_allclose(h5_dm.mz_axis, binary_dm.mz_axis, rtol=1e-7, atol=1e-10)
        np.testing.assert_allclose(zarr_dm.mz_axis, binary_dm.mz_axis, rtol=1e-7, atol=1e-10)

        for ion_idx in _sample_indices(binary_dm.n_ions):
            binary_image = binary_dm.get_ion_image(ion_idx)
            np.testing.assert_allclose(
                h5_dm.get_ion_image(ion_idx),
                binary_image,
                rtol=1e-5,
                atol=1e-8,
            )
            np.testing.assert_allclose(
                zarr_dm.get_ion_image(ion_idx),
                binary_image,
                rtol=1e-5,
                atol=1e-8,
            )
            np.testing.assert_allclose(
                h5_dm.get_ion_vector(ion_idx),
                binary_dm.get_ion_vector(ion_idx),
                rtol=1e-5,
                atol=1e-8,
            )
            np.testing.assert_allclose(
                zarr_dm.get_ion_vector(ion_idx),
                binary_dm.get_ion_vector(ion_idx),
                rtol=1e-5,
                atol=1e-8,
            )
