from typing import Optional
from massflow.preprocess.numba.numba_runtime import apply_numba_runtime
from massflow.tools import get_logger
from massflow.preprocess.numba.numba_runtime import get_global_numba_runtime, get_num_threads

logger = get_logger("massflow.tools.dm_process")


def flat_speed_process(
    flats,
    func,
    batch_kwargs: dict,
    numba_max_threads: Optional[int] = None,
) :

    apply_numba_runtime(override_workers=numba_max_threads)
    logger.info(f"numba treads set to {get_global_numba_runtime()['max_workers']}\r\n"
                f"numba system info: {get_num_threads()} threads")

    for batch_idx, batch in enumerate(flats, start=1):
        _ = func(flat=batch, **batch_kwargs)
