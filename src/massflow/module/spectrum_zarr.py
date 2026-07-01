"""Lazy spectrum backed by a MassFlow spectra Zarr store."""

from __future__ import annotations

from typing import Any

from massflow.module.spectrum import Spectrum


class SpectrumZarr(Spectrum):
    """Spectrum placeholder that reads mz/intensity arrays from Zarr on demand."""

    def __init__(
        self,
        reader: Any,
        index: int,
        coordinates,
        mz_list=None,
        shared_mz_list=None,
        sort_by_mz: bool = True,
    ):
        self._reader = reader
        self._index = int(index)
        super().__init__(
            mz_list=mz_list,
            intensity=None,
            coordinate=coordinates,
            shared_mz_list=shared_mz_list,
            sort_by_mz=sort_by_mz,
        )

    def resolve_data(self, **kwargs):
        """Load this spectrum from the owning Zarr reader."""
        if self._intensity is not None:
            return
        mz, intensity = self._reader.read_spectrum(self._index)
        if self.shared_mz_list is not None:
            self._mz_list = self.shared_mz_list
        else:
            self._mz_list = mz
        self._intensity = intensity
