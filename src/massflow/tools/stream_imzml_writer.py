"""
Main optimizations:
1. Remove SHA-1 hash calculation (use fixed placeholder)
2. Remove statistics calculation (min/max/base peak/TIC use 0 placeholders)
3. Optimize continuous mode (no m/z array comparison check)
4. Avoid unnecessary memory copying (zero-copy write strategy)
5. Remove rounding operation (rounding in compression mode)
"""

from __future__ import print_function

import os
import uuid as uuid_module
from collections import namedtuple, OrderedDict
from typing import Optional, List, Tuple, Union, Any
import numpy as np

from wheezy.template import Engine, CoreExtension, DictLoader
from pyimzml.compression import NoCompression

_Spectrum = namedtuple(
    '_Spectrum', 
    'coords mz_len mz_offset mz_enc_len int_len int_offset int_enc_len mz_min mz_max mz_base int_base int_tic userParams'
)

IMZML_TEMPLATE = """\
@require(uuid, sha1sum, mz_data_type, int_data_type, run_id, spectra, mode, obo_codes, obo_names, mz_compression, int_compression, polarity, spec_type, scan_direction, scan_pattern, scan_type, line_scan_direction)
<?xml version="1.0" encoding="ISO-8859-1"?>
<mzML xmlns="http://psi.hupo.org/ms/mzml" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://psi.hupo.org/ms/mzml http://psidev.info/files/ms/mzML/xsd/mzML1.1.0_idx.xsd" version="1.1">
  <cvList count="2">
    <cv uri="http://psidev.cvs.sourceforge.net/*checkout*/psidev/psi/psi-ms/mzML/controlledVocabulary/psi-ms.obo" fullName="Proteomics Standards Initiative Mass Spectrometry Ontology" id="MS" version="3.65.0"/>
    <cv uri="http://obo.cvs.sourceforge.net/*checkout*/obo/obo/ontology/phenotype/unit.obo" fullName="Unit Ontology" id="UO" version="12:10:2011"/>
  </cvList>

  <fileDescription>
    <fileContent>
      <cvParam cvRef="MS" accession="MS:1000579" name="MS1 spectrum" value=""/>
      @if spec_type=='centroid':
      <cvParam cvRef="MS" accession="MS:1000127" name="centroid spectrum" value=""/>
      @elif spec_type=='profile':
      <cvParam cvRef="MS" accession="MS:1000128" name="profile spectrum" value=""/>
      @end
      <cvParam cvRef="IMS" accession="IMS:@obo_codes[mode]" name="@mode" value=""/>
      <cvParam cvRef="IMS" accession="IMS:1000080" name="universally unique identifier" value="@uuid"/>
      <cvParam cvRef="IMS" accession="IMS:1000091" name="ibd SHA-1" value="@sha1sum"/>
    </fileContent>
  </fileDescription>

  <referenceableParamGroupList count="4">
    <referenceableParamGroup id="mzArray">
      <cvParam cvRef="MS" accession="MS:@obo_codes[mz_compression]" name="@mz_compression" value=""/>
      <cvParam cvRef="MS" accession="MS:1000514" name="m/z array" unitCvRef="MS" unitAccession="MS:1000040" unitName="m/z"/>
      <cvParam cvRef="MS" accession="MS:@obo_codes[mz_data_type]" name="@mz_data_type" value=""/>
      <cvParam cvRef="IMS" accession="IMS:1000101" name="external data" value="true"/>
    </referenceableParamGroup>
    <referenceableParamGroup id="intensityArray">
      <cvParam cvRef="MS" accession="MS:@obo_codes[int_data_type]" name="@int_data_type" value=""/>
      <cvParam cvRef="MS" accession="MS:1000515" name="intensity array" unitCvRef="MS" unitAccession="MS:1000131" unitName="number of detector counts"/>
      <cvParam cvRef="MS" accession="MS:@obo_codes[int_compression]" name="@int_compression" value=""/>
      <cvParam cvRef="IMS" accession="IMS:1000101" name="external data" value="true"/>
    </referenceableParamGroup>
    <referenceableParamGroup id="scan1">
      <cvParam cvRef="MS" accession="MS:1000093" name="increasing m/z scan"/>
      <cvParam cvRef="MS" accession="MS:1000512" name="filter string" value=""/>
    </referenceableParamGroup>
    <referenceableParamGroup id="spectrum1">
      <cvParam cvRef="MS" accession="MS:1000579" name="MS1 spectrum" value=""/>
      <cvParam cvRef="MS" accession="MS:1000511" name="ms level" value="0"/>
      @if spec_type=='centroid':
      <cvParam cvRef="MS" accession="MS:1000127" name="centroid spectrum" value=""/>
      @elif spec_type=='profile':
      <cvParam cvRef="MS" accession="MS:1000128" name="profile spectrum" value=""/>
      @end
      @if polarity=='positive':
      <cvParam cvRef="MS" accession="MS:1000130" name="positive scan" value=""/>
      @elif polarity=='negative':
      <cvParam cvRef="MS" accession="MS:1000129" name="negative scan" value=""/>
      @end
    </referenceableParamGroup>
  </referenceableParamGroupList>

  <softwareList count="1">
    <software id="pyimzml" version="0.0001">
      <cvParam cvRef="MS" accession="MS:1000799" name="custom unreleased software tool" value="pyimzml exporter"/>
    </software>
  </softwareList>

  <scanSettingsList count="1">
    <scanSettings id="scanSettings1">
      <cvParam cvRef="IMS" accession="IMS:@obo_codes[scan_direction]" name="@obo_names[scan_direction]" value=""/>
      <cvParam cvRef="IMS" accession="IMS:@obo_codes[scan_pattern]" name="@obo_names[scan_pattern]" value=""/>
      <cvParam cvRef="IMS" accession="IMS:@obo_codes[scan_type]" name="@obo_names[scan_type]" value=""/>
      <cvParam cvRef="IMS" accession="IMS:@obo_codes[line_scan_direction]" name="@obo_names[line_scan_direction]" value=""/>
      <cvParam cvRef="IMS" accession="IMS:1000042" name="max count of pixels x" value="@{(max(s.coords[0] for s in spectra))!!s}"/>
      <cvParam cvRef="IMS" accession="IMS:1000043" name="max count of pixels y" value="@{(max(s.coords[1] for s in spectra))!!s}"/>
    </scanSettings>
  </scanSettingsList>

  <instrumentConfigurationList count="1">
    <instrumentConfiguration id="IC1">
    </instrumentConfiguration>
  </instrumentConfigurationList>

  <dataProcessingList count="1">
    <dataProcessing id="export_from_pyimzml">
      <processingMethod order="0" softwareRef="pyimzml">
        <cvParam cvRef="MS" accession="MS:1000530" name="file format conversion" value="Output to imzML"/>
      </processingMethod>
    </dataProcessing>
  </dataProcessingList>

  <run defaultInstrumentConfigurationRef="IC1" id="@run_id">
    <spectrumList count="@{len(spectra)!!s}" defaultDataProcessingRef="export_from_pyimzml">
      @for index, s in enumerate(spectra):
      <spectrum defaultArrayLength="0" id="spectrum=@{(index+1)!!s}" index="@{(index+1)!!s}">
        <referenceableParamGroupRef ref="spectrum1"/>
        <cvParam cvRef="MS" accession="MS:1000528" name="lowest observed m/z" value="@{s.mz_min!!s}" unitCvRef="MS" unitAccession="MS:1000040" unitName="m/z"/>
        <cvParam cvRef="MS" accession="MS:1000527" name="highest observed m/z" value="@{s.mz_max!!s}" unitCvRef="MS" unitAccession="MS:1000040" unitName="m/z"/>
        <cvParam cvRef="MS" accession="MS:1000504" name="base peak m/z" value="@{s.mz_base!!s}" unitCvRef="MS" unitAccession="MS:1000040" unitName="m/z"/>
        <cvParam cvRef="MS" accession="MS:1000505" name="base peak intensity" value="@{s.int_base!!s}" unitCvRef="MS" unitAccession="MS:1000131" unitName="number of detector counts"/>
        <cvParam cvRef="MS" accession="MS:1000285" name="total ion current" value="@{s.int_tic!!s}"/>
        <scanList count="1">
          <cvParam cvRef="MS" accession="MS:1000795" name="no combination"/>
          <scan>
            <referenceableParamGroupRef ref="scan1"/>
            <cvParam accession="IMS:1000050" cvRef="IMS" name="position x" value="@{s.coords[0]!!s}"/>
            <cvParam accession="IMS:1000051" cvRef="IMS" name="position y" value="@{s.coords[1]!!s}"/>
            @if len(s.coords) == 3:
            <cvParam accession="IMS:1000052" cvRef="IMS" name="position z" value="@{s.coords[2]!!s}"/>
            @end
            @if s.userParams:
                @for up in s.userParams:
                <userParam name="@up['name']" value="@up['value']"/> 
                @end
            @end
          </scan>
        </scanList>
        <binaryDataArrayList count="2">
          <binaryDataArray encodedLength="0">
            <referenceableParamGroupRef ref="mzArray"/>
            <cvParam accession="IMS:1000103" cvRef="IMS" name="external array length" value="@{s.mz_len!!s}"/>
            <cvParam accession="IMS:1000104" cvRef="IMS" name="external encoded length" value="@{s.mz_enc_len!!s}"/>
            <cvParam accession="IMS:1000102" cvRef="IMS" name="external offset" value="@{s.mz_offset!!s}"/>
            <binary/>
          </binaryDataArray>
          <binaryDataArray encodedLength="0">
            <referenceableParamGroupRef ref="intensityArray"/>
            <cvParam accession="IMS:1000103" cvRef="IMS" name="external array length" value="@{s.int_len!!s}"/>
            <cvParam accession="IMS:1000104" cvRef="IMS" name="external encoded length" value="@{s.int_enc_len!!s}"/>
            <cvParam accession="IMS:1000102" cvRef="IMS" name="external offset" value="@{s.int_offset!!s}"/>
            <binary/>
          </binaryDataArray>
        </binaryDataArrayList>
      </spectrum>
      @end
    </spectrumList>
  </run>
</mzML>
"""


PLACEHOLDER_SHA1 = "0" * 40

class ImzMLWriter:
    """
    High-performance imzML writer - Drop-in Replacement for pyimzml.ImzMLWriter
    
    Extremely optimized for "Swap to Disk" scenarios, write speed approaching physical disk limits.
    
    Optimization strategies:
    1. Completely remove SHA-1/MD5 hash calculation (use fixed placeholder)
    2. Optionally skip statistics calculation (min/max/base peak/TIC)
    3. No m/z array comparison in continuous mode (trust caller)
    4. Zero-copy write strategy (avoid array copying during dtype conversion)
    5. Remove compression rounding operation
    """

    def __init__(
        self,
        output_filename: str,
        mz_dtype: np.dtype = np.float64, #type: ignore
        intensity_dtype: np.dtype = np.float32, #type: ignore
        mode: str = "auto",
        spec_type: str = "centroid",
        scan_direction: str = "top_down",
        line_scan_direction: str = "line_left_right",
        scan_pattern: str = "one_way",
        scan_type: str = "horizontal_line",
        mz_compression: Any = NoCompression(),
        intensity_compression: Any = NoCompression(),
        polarity: Optional[str] = None,
        compute_stats: bool = False,
    ):
        self.mz_dtype = np.dtype(mz_dtype)
        self.intensity_dtype = np.dtype(intensity_dtype)
        self.mode = mode
        self.spec_type = spec_type
        self.mz_compression = mz_compression
        self.intensity_compression = intensity_compression
        self.compute_stats = compute_stats

        # File path processing
        self.run_id = os.path.splitext(output_filename)[0]
        self.filename = self.run_id + ".imzML"
        self.ibd_filename = self.run_id + ".ibd"

        # Open files
        self.xml = open(self.filename, 'w') # pylint: disable=W1514
        self.ibd = open(self.ibd_filename, 'wb')

        # UUID
        self.uuid = uuid_module.uuid4()

        # Scan parameters
        self.scan_direction = scan_direction
        self.scan_pattern = scan_pattern
        self.scan_type = scan_type
        self.line_scan_direction = line_scan_direction

        # Write UUID to ibd file header
        self.ibd.write(self.uuid.bytes)

        # Initialize template engine
        self.wheezy_engine = Engine(
            loader=DictLoader({'imzml': IMZML_TEMPLATE}),
            extensions=[CoreExtension()]
        )
        self.imzml_template = self.wheezy_engine.get_template('imzml')

        # Spectrum data list
        self.spectra: List[_Spectrum] = []

        # Continuous mode: Store position information of the first m/z array
        self.first_mz: Optional[Tuple[int, int, int]] = None

        # LRU cache needed for auto mode (maintain compatibility)
        self._lru_cache: OrderedDict = OrderedDict()
        self._lru_maxlen = 10

        # Polarity setting
        self._setPolarity(polarity)

        # Pre-calculate dtype byte sizes
        self._mz_itemsize = self.mz_dtype.itemsize
        self._int_itemsize = self.intensity_dtype.itemsize

    @staticmethod
    def _np_type_to_name(dtype: np.dtype) -> str:
        """Convert numpy dtype to imzML specification name"""
        dtype = np.dtype(dtype)
        if dtype == np.float32:
            return "32-bit float"
        elif dtype == np.float64:
            return "64-bit float"
        elif dtype == np.int32:
            return "32-bit integer"
        elif dtype == np.int64:
            return "64-bit integer"
        else:
            raise ValueError(f"Unsupported dtype: {dtype}")

    def _setPolarity(self, polarity: Optional[str]) -> None: # pylint: disable=C0103
        """Set polarity"""
        if polarity is None or polarity.lower() in ('positive', 'negative'):
            self.polarity = polarity.lower() if polarity else None
        else:
            raise ValueError(f"Invalid polarity: {polarity}. Must be 'positive', 'negative', or None")

    def _write_binary_fast(self, data: np.ndarray, target_dtype: np.dtype) -> Tuple[int, int, int]:
        """
        Fast binary write - Core optimization method
        
        Optimization strategies:
        1. Use np.asarray instead of np.array to avoid unnecessary copying
        2. If dtype matches and array is C-contiguous, directly use data.data buffer
        3. Skip rounding operation in compression
        
        Returns: (offset, length, encoded_length)
        """
        offset = self.ibd.tell()

        # Smart type conversion: Create new array only when necessary
        if data.dtype != target_dtype:
            data = data.astype(target_dtype, copy=False)

        # Ensure C-contiguous (most numpy arrays already are)
        if not data.flags['C_CONTIGUOUS']:
            data = np.ascontiguousarray(data)

        # Get raw bytes
        raw_bytes = data.tobytes()

        # Compression (if needed)
        if not isinstance(self.mz_compression, NoCompression):
            raw_bytes = self.mz_compression.compress(raw_bytes)

        # Write directly without hash calculation
        self.ibd.write(raw_bytes)

        return offset, len(data), len(raw_bytes)

    def _write_mz(self, mzs: np.ndarray) -> Tuple[int, int, int]:
        """Write m/z array"""
        return self._write_binary_fast(mzs, self.mz_dtype)

    def _write_intensities(self, intensities: np.ndarray) -> Tuple[int, int, int]:
        """Write intensity array"""
        offset = self.ibd.tell()

        # Smart type conversion
        if intensities.dtype != self.intensity_dtype:
            intensities = intensities.astype(self.intensity_dtype, copy=False)

        if not intensities.flags['C_CONTIGUOUS']:
            intensities = np.ascontiguousarray(intensities)

        raw_bytes = intensities.tobytes()

        # Compression
        if not isinstance(self.intensity_compression, NoCompression):
            raw_bytes = self.intensity_compression.compress(raw_bytes)

        self.ibd.write(raw_bytes)

        return offset, len(intensities), len(raw_bytes)

    def addSpectrum( # pylint: disable=W0102, C0103
        self,
        mzs: np.ndarray,
        intensities: np.ndarray,
        coords: Union[Tuple[int, int], Tuple[int, int, int]],
        userParams: List = [] # pylint: disable=C0103
    ) -> None:
        """
        Add a mass spectrum
        
        Parameters:
            mzs: m/z array
            intensities: Intensity array
            coords: Coordinates (x, y) or (x, y, z)
            userParams: User-defined parameter list
        """
        # Ensure numpy arrays
        if not isinstance(mzs, np.ndarray):
            mzs = np.asarray(mzs, dtype=self.mz_dtype)
        if not isinstance(intensities, np.ndarray):
            intensities = np.asarray(intensities, dtype=self.intensity_dtype)

        # Process m/z array based on mode
        if self.mode == "continuous":
            if self.first_mz is None:
                # First spectrum: write m/z array
                self.first_mz = self._write_mz(mzs)
            mz_offset, mz_len, mz_enc_len = self.first_mz

        elif self.mode == "processed":
            # Each spectrum has independent m/z
            mz_offset, mz_len, mz_enc_len = self._write_mz(mzs)

        elif self.mode == "auto":
            # Auto mode: Use simplified cache strategy
            mz_key = (len(mzs), mzs[0] if len(mzs) > 0 else 0, mzs[-1] if len(mzs) > 0 else 0)
            if mz_key in self._lru_cache:
                mz_offset, mz_len, mz_enc_len = self._lru_cache[mz_key]
            else:
                mz_data = self._write_mz(mzs)
                self._lru_cache[mz_key] = mz_data
                # Simple LRU: Delete oldest when exceeding limit
                if len(self._lru_cache) > self._lru_maxlen:
                    self._lru_cache.popitem(last=False)
                mz_offset, mz_len, mz_enc_len = mz_data
        else:
            raise ValueError(f"Unknown mode: {self.mode}")

        # Write intensity array
        int_offset, int_len, int_enc_len = self._write_intensities(intensities)

        # Calculate or skip statistics
        if self.compute_stats:
            mz_min = float(np.min(mzs))
            mz_max = float(np.max(mzs))
            ix_max = np.argmax(intensities)
            mz_base = float(mzs[ix_max]) if len(mzs) > 0 else 0.0
            int_base = float(intensities[ix_max]) if len(intensities) > 0 else 0.0
            int_tic = float(np.sum(intensities))
        else:
            # Use placeholder values - greatly improves performance
            mz_min = 0.0
            mz_max = 0.0
            mz_base = 0.0
            int_base = 0.0
            int_tic = 0.0

        # Create spectrum record
        s = _Spectrum(
            coords=coords,
            mz_len=mz_len,
            mz_offset=mz_offset,
            mz_enc_len=mz_enc_len,
            int_len=int_len,
            int_offset=int_offset,
            int_enc_len=int_enc_len,
            mz_min=mz_min,
            mz_max=mz_max,
            mz_base=mz_base,
            int_base=int_base,
            int_tic=int_tic,
            userParams=userParams
        )
        self.spectra.append(s)

    def addMatrix( # pylint: disable=C0103
            self,
            mz_data: np.ndarray,
            intensity_matrix: np.ndarray,
            lengths: np.ndarray,
            coordinates: np.ndarray
        ) -> None:
        """
        swap matrix data -> disk
        """
        if intensity_matrix.shape[0] == 0:
            return

        intensity_matrix = self._ensure_matrix(intensity_matrix, self.intensity_dtype)
        mz_data = self._ensure_matrix(mz_data, self.mz_dtype)

        lengths_arr = np.asarray(lengths, dtype=np.int64)
        coord_fn = self._make_coord_fn(coordinates)

        is_continuous = (self.mode == "continuous") or (mz_data.ndim == 1)
        self.first_mz = self._write_mz(mz_data) if is_continuous and self.first_mz is None else self.first_mz
        shared_mz = self.first_mz if self.first_mz is not None else (0, 0, 0)

        if self._can_use_fast_path(is_continuous, lengths_arr, intensity_matrix.shape[1]):
            self._write_matrix_fast(
                intensity_matrix=intensity_matrix,
                batch_size=int(intensity_matrix.shape[0]),
                max_length=int(intensity_matrix.shape[1]),
                coord_fn=coord_fn,
                shared_mz=shared_mz
            )
            return

        self._write_matrix_rows(
            mz_data=mz_data,
            intensity_matrix=intensity_matrix,
            lengths_arr=lengths_arr,
            batch_size=int(intensity_matrix.shape[0]),
            coord_fn=coord_fn,
            is_continuous=is_continuous,
            shared_mz=shared_mz
        )
    
    def _ensure_matrix(self, arr: np.ndarray, dtype: np.dtype) -> np.ndarray:
        if arr.dtype != dtype:
            arr = arr.astype(dtype, copy=False)
        if not arr.flags["C_CONTIGUOUS"]:
            arr = np.ascontiguousarray(arr)
        return arr

    def _make_coord_fn(self, coordinates: np.ndarray):
        if coordinates.shape[1] == 2:
            def coord_tuple(i: int):
                return int(coordinates[i, 0]), int(coordinates[i, 1])
        else:
            def coord_tuple(i: int):
                return int(coordinates[i, 0]), int(coordinates[i, 1]), int(coordinates[i, 2])
        return coord_tuple

    def _can_use_fast_path(self, is_continuous: bool, lengths_arr: np.ndarray, max_length: int) -> bool:
        if not is_continuous or int(lengths_arr[0]) != max_length:
            return False
        if not isinstance(self.intensity_compression, NoCompression):
            return False
        return bool(np.all(lengths_arr == lengths_arr[0]))

    def _write_matrix_fast(
            self,
            intensity_matrix: np.ndarray,
            batch_size: int,
            max_length: int,
            coord_fn,
            shared_mz: tuple[int, int, int],
        ) -> None:
        shared_mz_offset, shared_mz_len, shared_mz_enc_len = shared_mz

        row_bytes = max_length * self._int_itemsize
        base_offset = self.ibd.tell()

        self.ibd.write(memoryview(intensity_matrix))

        offsets = base_offset + np.arange(batch_size, dtype=np.int64) * row_bytes

        self.spectra.extend([
            self._get_spectrum(
                coords=coord_fn(i),
                mz_len=shared_mz_len,
                mz_offset=shared_mz_offset,
                mz_enc_len=shared_mz_enc_len,
                int_len=max_length,
                int_offset=int(offsets[i]),
                int_enc_len=row_bytes,
            )
            for i in range(batch_size)
        ])

    def _write_matrix_rows(
            self,
            mz_data: np.ndarray,
            intensity_matrix: np.ndarray,
            lengths_arr: np.ndarray,
            batch_size: int,
            coord_fn,
            is_continuous: bool,
            shared_mz: tuple[int, int, int],
        ) -> None:
        no_int_compress = isinstance(self.intensity_compression, NoCompression)
        no_mz_compress = isinstance(self.mz_compression, NoCompression)
        use_mz_cache = (self.mode == "auto" and not is_continuous)

        shared_mz_offset, shared_mz_len, shared_mz_enc_len = shared_mz

        spectra = []
        for i in range(batch_size):
            length = int(lengths_arr[i])
            coords = coord_fn(i)

            if length <= 0:
                cur_off = self.ibd.tell()
                spectra.append(self._get_spectrum(coords=coords, mz_offset=cur_off, int_offset=cur_off))
                continue

            if is_continuous:
                mz_off, mz_len, mz_enc = shared_mz_offset, shared_mz_len, shared_mz_enc_len
            else:
                mz_off, mz_len, mz_enc = self._write_mz_row(
                    mz_data=mz_data,
                    i=i,
                    length=length,
                    use_cache=use_mz_cache,
                    no_mz_compress=no_mz_compress,
                )

            int_off, int_enc_len = self._write_int_row(
                intensity_matrix=intensity_matrix,
                i=i,
                length=length,
                no_int_compress=no_int_compress,
            )

            spectra.append(self._get_spectrum(
                coords=coords,
                mz_len=mz_len,
                mz_offset=mz_off,
                mz_enc_len=mz_enc,
                int_len=length,
                int_offset=int_off,
                int_enc_len=int_enc_len,
            ))

        self.spectra.extend(spectra)

    def _mz_cache_key(self, mz_row: np.ndarray, length: int):
        return (
            length,
            float(mz_row[0]) if length > 0 else 0,
            float(mz_row[length // 2]) if length >= 3 else 0,
            float(mz_row[-1]) if length > 0 else 0,
            float(np.sum(mz_row)),
        )

    def _write_mz_row(self, mz_data: np.ndarray, i: int, length: int, use_cache: bool, no_mz_compress: bool):
        mz_row = mz_data[i, :length]

        if use_cache:
            key = self._mz_cache_key(mz_row, length)
            hit = self._lru_cache.get(key)
            if hit is not None:
                return hit

        mz_off = self.ibd.tell()
        mz_buf = mz_row.tobytes()
        if not no_mz_compress:
            mz_buf = self.mz_compression.compress(mz_buf)
        self.ibd.write(mz_buf)
        mz_len = length
        mz_enc = len(mz_buf)

        if use_cache:
            self._lru_cache[key] = (mz_off, mz_len, mz_enc)
            if len(self._lru_cache) > self._lru_maxlen:
                self._lru_cache.popitem(last=False)

        return (mz_off, mz_len, mz_enc)

    def _write_int_row(self, intensity_matrix: np.ndarray, i: int, length: int, no_int_compress: bool):
        int_off = self.ibd.tell()
        int_buf = intensity_matrix[i, :length].tobytes()
        if not no_int_compress:
            int_buf = self.intensity_compression.compress(int_buf)
        self.ibd.write(int_buf)
        return int_off, len(int_buf)

    @staticmethod
    def _get_spectrum(coords=0, mz_len=0, mz_offset=0, mz_enc_len=0, int_len=0, int_offset=0, int_enc_len=0) -> _Spectrum:
        """
        get _Spectrum
        """
        return _Spectrum(
            coords=coords,
            mz_len=mz_len,
            mz_offset=mz_offset,
            mz_enc_len=mz_enc_len,
            int_len=int_len,
            int_offset=int_offset,
            int_enc_len=int_enc_len,
            mz_min=0.0,
            mz_max=0.0,
            mz_base=0.0,
            int_base=0.0,
            int_tic=0.0,
            userParams=[],
        )

    def _write_xml(self) -> None:
        """Write imzML XML file"""
        # The following variables are passed to the template engine via locals(), Pylint cannot recognize
        # pylint: disable=W0641
        spectra = self.spectra
        mz_data_type = self._np_type_to_name(self.mz_dtype)
        int_data_type = self._np_type_to_name(self.intensity_dtype)

        obo_codes = {
            "32-bit integer": "1000519", 
            "16-bit float": "1000520",
            "32-bit float": "1000521",
            "64-bit integer": "1000522",
            "64-bit float": "1000523",
            "continuous": "1000030",
            "processed": "1000031",
            "zlib compression": "1000574",
            "no compression": "1000576",
            "line_bottom_up": "1000492",
            "line_left_right": "1000491",
            "line_right_left": "1000490",
            "line_top_down": "1000493",
            "bottom_up": "1000400",
            "left_right": "1000402",
            "right_left": "1000403",
            "top_down": "1000401",
            "meandering": "1000410",
            "one_way": "1000411",
            "random_access": "1000412",
            "horizontal_line": "1000480",
            "vertical_line": "1000481"
        }
        obo_names = {
            "line_bottom_up": "linescan bottom up",
            "line_left_right": "linescan left right",
            "line_right_left": "linescan right left",
            "line_top_down": "linescan top down",
            "bottom_up": "bottom up",
            "left_right": "left right",
            "right_left": "right left",
            "top_down": "top down",
            "meandering": "meandering",
            "one_way": "one way",
            "random_access": "random access",
            "horizontal_line": "horizontal line scan",
            "vertical_line": "vertical line scan"
        }

        uuid = ("{%s}" % self.uuid).upper()
        sha1sum = PLACEHOLDER_SHA1  # 使用固定占位符，跳过计算
        run_id = self.run_id

        # 确定模式
        if self.mode == 'auto':
            mode = "processed" if len(self._lru_cache) > 1 else "continuous"
        else:
            mode = self.mode

        spec_type = self.spec_type
        mz_compression = self.mz_compression.name
        int_compression = self.intensity_compression.name
        polarity = self.polarity
        scan_direction = self.scan_direction
        scan_pattern = self.scan_pattern
        scan_type = self.scan_type
        line_scan_direction = self.line_scan_direction

        self.xml.write(self.imzml_template.render(locals()))

    def close(self) -> None:
        """
        关闭写入器，写入 XML 文件并关闭所有文件句柄。
        使用 with 语句时会自动调用。
        """
        self.finish()

    def finish(self) -> None:
        """close() 的别名"""
        self.ibd.close()
        self._write_xml()
        self.xml.close()

    def __enter__(self) -> 'ImzMLWriter':
        return self

    def __exit__(self, exc_t, exc_v, trace) -> None:
        if exc_t is None:
            self.finish()
        else:
            self.ibd.close()
            self.xml.close()
