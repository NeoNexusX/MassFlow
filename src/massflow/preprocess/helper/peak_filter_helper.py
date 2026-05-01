from typing import Optional

import numpy as np

# 总接口定义
def peak_filter(
    mz_data: np.ndarray,
    intensity: np.ndarray,
    lengths: np.ndarray,
    *,
    method: str,
    counts: np.ndarray,
    min_count: Optional[int] = None,
    min_frequency: Optional[float] = None,
):
    _ = mz_data, intensity, lengths, method, counts, min_count, min_frequency
