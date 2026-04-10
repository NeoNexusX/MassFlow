"""Utilities for MassFlow.

Note: This package coexists with a legacy `__init___.py` file in the same
folder. Use `massflow.tools` imports going forward.
"""

from .logger import get_logger  # re-export for convenience
from .funs import is_valid_file, dispatch_with_supported_kwargs,prepare_flat_inputs, lengths_to_offsets  # re-export for convenience
