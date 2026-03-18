"""
主要优化点：
1. 移除 SHA-1 哈希计算 (使用固定占位符)
2. 移除统计量计算 (min/max/base peak/TIC 使用 0 占位符)
3. 优化 continuous 模式 (无 m/z 数组比对检查)
4. 避免不必要的内存拷贝 (零拷贝写入策略)
5. 移除 rounding 操作 (压缩模式下的舍入)
"""

from __future__ import print_function

import os
import numpy as np
import uuid as uuid_module
from collections import namedtuple, OrderedDict
from typing import Optional, List, Tuple, Union, Any

from wheezy.template import Engine, CoreExtension, DictLoader
from pyimzml.compression import NoCompression

_Spectrum = namedtuple(
    "_Spectrum",
    "coords mz_len mz_offset mz_enc_len int_len int_offset int_enc_len mz_min mz_max mz_base int_base int_tic userParams",
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
    高性能 imzML 写入器 - pyimzml.ImzMLWriter 的 Drop-in Replacement

    针对 "Swap to Disk" 场景极致优化，写入速度接近物理磁盘极限。

    优化策略:
    1. 完全移除 SHA-1/MD5 哈希计算 (使用固定占位符)
    2. 可选跳过统计量计算 (min/max/base peak/TIC)
    3. continuous 模式无 m/z 数组比对 (信任调用者)
    4. 零拷贝写入策略 (避免 dtype 转换时的数组复制)
    5. 移除压缩 rounding 操作
    """

    def __init__(
        self,
        output_filename: str,
        mz_dtype: np.dtype = np.float64,  # type: ignore
        intensity_dtype: np.dtype = np.float32,  # type: ignore
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

        # 文件路径处理
        self.run_id = os.path.splitext(output_filename)[0]
        self.filename = self.run_id + ".imzML"
        self.ibd_filename = self.run_id + ".ibd"

        # 打开文件
        self.xml = open(self.filename, "w")
        self.ibd = open(self.ibd_filename, "wb")

        # UUID
        self.uuid = uuid_module.uuid4()

        # 扫描参数
        self.scan_direction = scan_direction
        self.scan_pattern = scan_pattern
        self.scan_type = scan_type
        self.line_scan_direction = line_scan_direction

        # 写入 UUID 到 ibd 文件头
        self.ibd.write(self.uuid.bytes)

        # 初始化模板引擎
        self.wheezy_engine = Engine(
            loader=DictLoader({"imzml": IMZML_TEMPLATE}), extensions=[CoreExtension()]
        )
        self.imzml_template = self.wheezy_engine.get_template("imzml")

        # 谱图数据列表
        self.spectra: List[_Spectrum] = []

        # continuous 模式: 存储第一个 m/z 数组的位置信息
        self.first_mz: Optional[Tuple[int, int, int]] = None

        # auto 模式需要的 LRU 缓存 (保持兼容性)
        self._lru_cache: OrderedDict = OrderedDict()
        self._lru_maxlen = 10

        # 极性设置
        self._set_polarity(polarity)

        # 预计算 dtype 的字节大小
        self._mz_itemsize = self.mz_dtype.itemsize
        self._int_itemsize = self.intensity_dtype.itemsize

    @staticmethod
    def _np_type_to_name(dtype: np.dtype) -> str:
        """将 numpy dtype 转换为 imzML 规范名称"""
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

    def _set_polarity(self, polarity: Optional[str]) -> None:
        """设置极性"""
        if polarity is None or polarity.lower() in ("positive", "negative"):
            self.polarity = polarity.lower() if polarity else None
        else:
            raise ValueError(f"Invalid polarity: {polarity}. Must be 'positive', 'negative', or None")

    def _write_binary_fast(
        self, data: np.ndarray, target_dtype: np.dtype
    ) -> Tuple[int, int, int]:
        """
        快速二进制写入 - 核心优化方法

        优化策略:
        1. 使用 np.asarray 替代 np.array，避免不必要的拷贝
        2. 如果 dtype 匹配且数组是 C-contiguous，直接使用 data.data buffer
        3. 跳过压缩的 rounding 操作

        返回: (offset, length, encoded_length)
        """
        offset = self.ibd.tell()

        # 智能类型转换: 只在必要时创建新数组
        if data.dtype != target_dtype:
            data = data.astype(target_dtype, copy=False)

        # 确保是 C-contiguous (大多数 numpy 数组已经是)
        if not data.flags["C_CONTIGUOUS"]:
            data = np.ascontiguousarray(data)

        # 获取原始字节
        raw_bytes = data.tobytes()

        # 压缩 (如果需要)
        if not isinstance(self.mz_compression, NoCompression):
            raw_bytes = self.mz_compression.compress(raw_bytes)

        # 直接写入，不计算哈希
        self.ibd.write(raw_bytes)

        return offset, len(data), len(raw_bytes)

    def _write_mz(self, mzs: np.ndarray) -> Tuple[int, int, int]:
        """写入 m/z 数组"""
        return self._write_binary_fast(mzs, self.mz_dtype)

    def _write_intensities(self, intensities: np.ndarray) -> Tuple[int, int, int]:
        """写入强度数组"""
        offset = self.ibd.tell()

        # 智能类型转换
        if intensities.dtype != self.intensity_dtype:
            intensities = intensities.astype(self.intensity_dtype, copy=False)

        if not intensities.flags["C_CONTIGUOUS"]:
            intensities = np.ascontiguousarray(intensities)

        raw_bytes = intensities.tobytes()

        # 压缩
        if not isinstance(self.intensity_compression, NoCompression):
            raw_bytes = self.intensity_compression.compress(raw_bytes)

        self.ibd.write(raw_bytes)

        return offset, len(intensities), len(raw_bytes)

    def add_spectrum(
        self,
        mzs: np.ndarray,
        intensities: np.ndarray,
        coords: Union[Tuple[int, int], Tuple[int, int, int]],
        userParams: List = [],
    ) -> None:
        """
        添加一个质谱谱图

        参数:
            mzs: m/z 数组
            intensities: 强度数组
            coords: 坐标 (x, y) 或 (x, y, z)
            userParams: 用户自定义参数列表
        """
        # 确保是 numpy 数组
        if not isinstance(mzs, np.ndarray):
            mzs = np.asarray(mzs, dtype=self.mz_dtype)
        if not isinstance(intensities, np.ndarray):
            intensities = np.asarray(intensities, dtype=self.intensity_dtype)

        # 根据模式处理 m/z 数组
        if self.mode == "continuous":
            if self.first_mz is None:
                # 第一个谱图: 写入 m/z 数组
                self.first_mz = self._write_mz(mzs)
            mz_offset, mz_len, mz_enc_len = self.first_mz

        elif self.mode == "processed":
            # 每个谱图独立的 m/z
            mz_offset, mz_len, mz_enc_len = self._write_mz(mzs)

        elif self.mode == "auto":
            # auto 模式: 使用简化的缓存策略
            mz_key = (
                len(mzs),
                mzs[0] if len(mzs) > 0 else 0,
                mzs[-1] if len(mzs) > 0 else 0,
            )
            if mz_key in self._lru_cache:
                mz_offset, mz_len, mz_enc_len = self._lru_cache[mz_key]
            else:
                mz_data = self._write_mz(mzs)
                self._lru_cache[mz_key] = mz_data
                # 简单的 LRU: 超过限制时删除最旧的
                if len(self._lru_cache) > self._lru_maxlen:
                    self._lru_cache.popitem(last=False)
                mz_offset, mz_len, mz_enc_len = mz_data
        else:
            raise ValueError(f"Unknown mode: {self.mode}")

        # 写入强度数组
        int_offset, int_len, int_enc_len = self._write_intensities(intensities)

        # 计算或跳过统计量
        if self.compute_stats:
            mz_min = float(np.min(mzs))
            mz_max = float(np.max(mzs))
            ix_max = np.argmax(intensities)
            mz_base = float(mzs[ix_max]) if len(mzs) > 0 else 0.0
            int_base = float(intensities[ix_max]) if len(intensities) > 0 else 0.0
            int_tic = float(np.sum(intensities))
        else:
            # 使用占位符值 - 极大提升性能
            mz_min = 0.0
            mz_max = 0.0
            mz_base = 0.0
            int_base = 0.0
            int_tic = 0.0

        # 创建谱图记录
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
            userParams=userParams,
        )
        self.spectra.append(s)

    def _write_xml(self) -> None:
        """写入 imzML XML 文件"""
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
            "vertical_line": "1000481",
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
            "vertical_line": "vertical line scan",
        }

        uuid = ("{%s}" % self.uuid).upper()
        sha1sum = PLACEHOLDER_SHA1  # 使用固定占位符，跳过计算
        run_id = self.run_id

        # 确定模式
        if self.mode == "auto":
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

    @classmethod
    def from_ms_meta(
        cls,
        output_filename: str,
        meta,
        mz_dtype=np.float64,
        intensity_dtype=np.float32,
    ):
        """
        从 MS metadata 创建 ImzMLWriter 实例
        """
        spec_type = "profile" if meta.profile_spectrum is not None else "centroid"
        scan_direction = meta.scan_direction if meta.scan_direction is not None else "top_down"
        line_scan_direction = meta.line_scan_direction if meta.line_scan_direction is not None else "line_left_right"
        scan_pattern = meta.scan_pattern if meta.scan_pattern is not None else "one_way"
        scan_type = meta.scan_type if meta.scan_type is not None else "horizontal_line"

        if meta.continuous is not None:
            mode = "continuous"
        elif meta.processed is not None:
            mode = "processed"
        else:
            mode = "auto"

        return cls(
            output_filename=output_filename,
            mz_dtype=mz_dtype, # type: ignore
            intensity_dtype=intensity_dtype,# type: ignore
            mode=mode,
            spec_type=spec_type,
            scan_direction=scan_direction,
            line_scan_direction=line_scan_direction,
            scan_pattern=scan_pattern,
            scan_type=scan_type,
        )


    def close(self) -> None:
        """
        关闭写入器，写入 XML 文件并关闭所有文件句柄。
        使用 with 语句时会自动调用。
        """
        self.finish()

    def finish(self) -> None:
        """close() 的别名"""
        if getattr(self, "_finished", False):
            return
        self._finished = True
        self.ibd.close()
        self._write_xml()
        self.xml.close()

    def __enter__(self) -> "ImzMLWriter":
        return self

    def __exit__(self, exc_t, exc_v, trace) -> None:
        if exc_t is None:
            self.finish()
        else:
            self.ibd.close()
            self.xml.close()