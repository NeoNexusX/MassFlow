"""
Tests for filter_helper module, specifically for the smooth_signal_loess function.
"""
import pytest
import numpy as np
from massflow.preprocess.filter_helper import (
    smooth_signal_loess,
    smoother
)
from massflow.preprocess.peak_align_helper import _smooth_loess_tricubic


class TestSmoothSignalLoess:
    """Tests for the smooth_signal_loess function."""

    def test_basic_smoothing(self):
        """Test basic smoothing functionality."""
        y = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 
                      11.0, 12.0, 13.0, 14.0, 15.0])
        result = smooth_signal_loess(y, window=5)
        assert result.shape == y.shape
        assert isinstance(result, np.ndarray)

    def test_backward_compatibility(self):
        """Test that _smooth_loess_tricubic is a backward-compatible alias."""
        y = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0])
        result1 = smooth_signal_loess(y, window=5)
        result2 = _smooth_loess_tricubic(y, window=5)
        np.testing.assert_array_almost_equal(result1, result2)

    def test_smoother_loess_method(self):
        """Test that smoother function supports loess method."""
        y = np.array([1.0, 2.0, 1.5, 3.0, 2.5, 4.0, 3.5, 5.0, 4.5, 6.0])
        result = smoother(y, method='loess', window=5)
        assert result.shape == y.shape
        assert isinstance(result, np.ndarray)

    def test_invalid_input_empty_array(self):
        """Test that empty array raises ValueError."""
        with pytest.raises(ValueError, match="must be a non-empty 1D array"):
            smooth_signal_loess(np.array([]))

    def test_invalid_input_2d_array(self):
        """Test that 2D array raises ValueError."""
        with pytest.raises(ValueError, match="must be a non-empty 1D array"):
            smooth_signal_loess(np.array([[1, 2], [3, 4]]))

    def test_invalid_window_negative(self):
        """Test that negative window raises ValueError."""
        y = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
        with pytest.raises(ValueError, match="window must be a positive integer"):
            smooth_signal_loess(y, window=-1)

    def test_invalid_window_zero(self):
        """Test that zero window raises ValueError."""
        y = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
        with pytest.raises(ValueError, match="window must be a positive integer"):
            smooth_signal_loess(y, window=0)

    def test_window_larger_than_data(self):
        """Test handling when window is larger than data length."""
        y = np.array([1.0, 2.0, 3.0])
        result = smooth_signal_loess(y, window=100)
        np.testing.assert_array_equal(result, y)

    def test_default_window(self):
        """Test that default window of 10 works correctly."""
        y = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0,
                      11.0, 12.0, 13.0, 14.0, 15.0])
        result = smooth_signal_loess(y)
        assert result.shape == y.shape
