from typing import Any, Callable, Iterable

from massflow.concurrency import apply_compute_mask, phase
from massflow.tools import get_logger

logger = get_logger("massflow.tools.dm_process")


def flat_speed_process(
    flats: Iterable[Any],
    func: Callable[..., Any],
    batch_kwargs: dict[str, Any],
) -> None:
    """Process flat batches within the global compute budget."""
    with phase("compute"):
        applied = apply_compute_mask()
        logger.info("numba threads set to %s", applied)
        for batch in flats:
            func(flat=batch, **batch_kwargs)
