import inspect

from massflow.preprocess.batch_pre_fun import BatchPreprocess
from massflow.preprocess.preprocessor import Preprocessor


def _without_receiver(params: list[inspect.Parameter], receiver_name: str) -> list[inspect.Parameter]:
    return [p for p in params if p.name != receiver_name]


def _normalize_preprocessor_params(params: list[inspect.Parameter]) -> list[inspect.Parameter]:
    return _without_receiver(params, "self")


def _normalize_batch_params(params: list[inspect.Parameter]) -> list[inspect.Parameter]:
    return _without_receiver(params, "batch_spectra")


def _assert_signature_compatible(pre_fn, batch_fn) -> None:
    pre_params = _normalize_preprocessor_params(list(inspect.signature(pre_fn).parameters.values()))
    batch_params = _normalize_batch_params(list(inspect.signature(batch_fn).parameters.values()))

    assert [p.name for p in pre_params] == [p.name for p in batch_params]
    assert [p.default for p in pre_params] == [p.default for p in batch_params]


def test_preprocessor_api_signatures_match_batch_preprocess() -> None:
    pairs = [
        (Preprocessor.baseline_correction, BatchPreprocess.baseline_correction_batch),
        (Preprocessor.noise_reduction, BatchPreprocess.noise_reduction_batch),
        (Preprocessor.normalization, BatchPreprocess.normalization_batch),
        (Preprocessor.peak_align, BatchPreprocess.peak_align_batch),
        (Preprocessor.peak_pick, BatchPreprocess.peak_pick_batch),
    ]

    for pre_fn, batch_fn in pairs:
        _assert_signature_compatible(pre_fn, batch_fn)
