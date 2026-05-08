"""
Preprocessor
  |
  |-- TaskPlanner
  |     |-- filter unsupported tasks by spectrum type
  |     |-- group consecutive batch tasks
  |     |-- mark final output segment
  |
  |-- PipelineRunner
        |
        |-- FlatPipelineExecutor
        |
        |-- StageFactory
        |     |-- FlatReaderStage - stage-a
        |     |-- BatchTaskStage - stage-b
        |     |-- FlatWriterStage - stage-c
        |
        |-- DatasetTaskExecutor
        |
        |-- OutputSinkFactory
              |-- ImzMLSink
              |-- IonImageSink
"""

from __future__ import annotations

from typing import Any, Self

from massflow.data_manager import MSDataManager
from massflow.preprocess.api import PreprocessorAPI, TaskScope
from massflow.tools.infer_spectrum_type import resolve_spectrum_type

from massflow.preprocess.pipeline.models import OutputFormat, PipelineConfig, PreprocessTask, TaskCallable
from massflow.preprocess.pipeline.planner import TaskPlanner
from massflow.preprocess.pipeline.runner import PipelineRunner


class PipelinePreprocessor(PreprocessorAPI):
    """Preprocessor implementation backed by the modular pipeline package."""

    def __init__(# pylint: disable=super-init-not-called
        self,
        data_manager: MSDataManager,
        *,
        batch_size: int = 256,
        queue_ab_size: int = 1,
        queue_bc_size: int = 1,
        temp_dir: str | None = "./temp",
        numba_max_threads: int | None = None,
        output_format: OutputFormat = "imzml",
        output_path: str | None = None,
        output_max_workers: int = 2,
    ):
        if data_manager is None:
            raise ValueError("data_manager must be provided.")

        self.data_manager = data_manager
        self.config = PipelineConfig(
            batch_size=batch_size,
            queue_ab_size=queue_ab_size,
            queue_bc_size=queue_bc_size,
            temp_dir=temp_dir,
            numba_max_threads=numba_max_threads,
            output_format=output_format,
            output_path=output_path,
            output_max_workers=output_max_workers,
        )
        self.config.validate()
        self._tasks: list[PreprocessTask] = []

    def _register_task(
        self,
        name: str,
        *,
        scope: TaskScope = "batch",
        apply_fn: TaskCallable,
        requires_dataset_prepare: bool = False,
        **kwargs: Any,
    ) -> Self:
        self._tasks.append(
            PreprocessTask(
                name=name,
                scope=scope,
                apply_fn=apply_fn,
                requires_dataset_prepare=requires_dataset_prepare,
                kwargs=kwargs,
            )
        )
        return self

    def start(self) -> Any:
        input_spectrum_type = resolve_spectrum_type(self.data_manager)
        plan = TaskPlanner(
            self._tasks,
            output_format=self.config.output_format,
        ).plan(input_spectrum_type)
        return PipelineRunner(self.config).run(data_manager=self.data_manager, plan=plan)
