from __future__ import annotations

from typing import Any, Callable, Optional

from massflow.data_manager import (
    MSDataManagerImzML,
    MSDataManagerZarr,
    IonImageDataManagerBinary,
    IonImageDataManagerH5,
    IonImageDataManagerZarr,
)
from massflow.preprocess.api import PreprocessorAPI
from massflow.preprocess.flat_pre_fun import FlatBatchResult
from massflow.preprocess.pipeline import (
    ION_IMAGE_OUTPUT_FORMATS,
    OUTPUT_FORMAT_IMZML,
    OUTPUT_FORMATS,
    InputMetaData,
    OutputFormat,
    PreprocessTask,
    TaskScope,
    imzml_to_ion_image_data_manager,
    make_task_segments,
    run_segments,
)
from massflow.tools.infer_spectrum_type import resolve_spectrum_type
from massflow.tools.logger import get_logger

logger = get_logger("massflow.preprocess.preprocessor")



class Preprocessor(PreprocessorAPI):
    """
    Hybrid preprocessing executor.

    - batch tasks are executed by 3-stage async pipeline
    - dataset tasks are executed on whole data manager (for peak align and cardinal backends)
    """

    def __init__(
        self,
        data_manager: MSDataManagerImzML,
        *,
        batch_size: int = 256,
        ion_block: int = 256,
        queue_ab_size: int = 1,
        queue_bc_size: int = 1,
        temp_dir: str = "./temp",
        numba_max_threads: Optional[int] = None,
        output_format: OutputFormat = OUTPUT_FORMAT_IMZML,
        output_path: str | None = None,
        output_max_workers: int = 2,
        compute_mean_spectrum: bool = False,
        compute_tic: bool = False,
    ): # pylint: disable=super-init-not-called
        if data_manager is None or batch_size <= 0 or batch_size > 9056 or queue_ab_size <= 0 or queue_bc_size <= 0 :
            logger.error(f"Invalid parameter values. please check:"
                         f"data_manager={data_manager}, batch_size={batch_size},"
                         f"queuesize={queue_ab_size} & {queue_bc_size}? <0?")
            raise ValueError("Invalid parameter values.")

        if output_format not in OUTPUT_FORMATS:
            raise ValueError(
                "output_format must be 'imzml', 'spectra_zarr', "
                "'ion_image_binary', 'ion_image_h5', or 'ion_image_zarr'."
            )
        self.data_manager = data_manager
        self.batch_size = batch_size
        self.temp_dir = temp_dir
        self.queue_ab_size = queue_ab_size
        self.queue_bc_size = queue_bc_size
        self.numba_max_threads = numba_max_threads
        self.output_format = output_format
        self.output_path = output_path
        self.output_max_workers = output_max_workers
        self.ion_block = ion_block
        self.compute_mean_spectrum = compute_mean_spectrum
        self.compute_tic = compute_tic

        self._tasks: list[PreprocessTask] = []

    def _register_task(
        self,
        name: str,
        *,
        scope: TaskScope = "batch",
        apply_fn: Callable[..., FlatBatchResult] | Callable[..., Any],
        **kwargs: Any,
    ) -> "Preprocessor":
        self._tasks.append(
            PreprocessTask(
                name=name,
                scope=scope,
                apply_fn=apply_fn,
                kwargs=kwargs,
            )
        )
        return self

    def start(self) -> MSDataManagerImzML | MSDataManagerZarr | IonImageDataManagerBinary | IonImageDataManagerH5 | IonImageDataManagerZarr:
        task_segments, input_meta = make_task_segments(
            list(self._tasks),
            input_spectrum_type=resolve_spectrum_type(self.data_manager),
            input_meta=InputMetaData(
                batch_size=self.batch_size,
                temp_dir=self.temp_dir,
                numba_max_threads=self.numba_max_threads,
                output_max_workers=self.output_max_workers,
                ion_block=self.ion_block,
                queue_ab_size=self.queue_ab_size,
                queue_bc_size=self.queue_bc_size,
                result_path=self.output_path,
                result_format=self.output_format, # type: ignore
                compute_mean_spectrum=self.compute_mean_spectrum,
                compute_tic=self.compute_tic,
            ),
        )
        if not task_segments and self.output_format in ION_IMAGE_OUTPUT_FORMATS:
            logger.warning("No preprocessing tasks registered, output will be written without processing.")
            return imzml_to_ion_image_data_manager(
                imzml_manager=self.data_manager, # type: ignore
                input_meta=input_meta,
                output_format=self.output_format,
                file_path=input_meta.result_path # type: ignore
            )

        return run_segments(self.data_manager, task_segments)
