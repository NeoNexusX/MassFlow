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
