from contextlib import contextmanager
import pyimzml.ImzMLParser
from massflow.tools.logger import get_logger

logger = get_logger("imzml_patch")

@contextmanager
def risky_imzml_loader():
    """
    上下文管理器：临时给 pyimzml 打补丁，使其能容忍缺失的 scanSettings 元数据。

    使用方法:
        with risky_imzml_loader():
            data_manager.load_full_data_from_file()

    原理:
        pyimzml 的 _get_cv_param 函数在遇到 None 元素时会报错 (AttributeError)。
        在进入上下文时替换该函数，遇到 None 直接返回 None，从而绕过报错。
        退出上下文时自动还原。
    """

    # 1. 保存原版函数引用
    original_get_cv_param = pyimzml.ImzMLParser._get_cv_param # pylint: disable=protected-access

    # 2. 定义临时函数
    def tolerant_get_cv_param(elem, accession, deep=False, convert=False):
        if elem is None:
            return None
        return original_get_cv_param(elem, accession, deep, convert)

    try:
        pyimzml.ImzMLParser._get_cv_param = tolerant_get_cv_param # pylint: disable=protected-access
        yield
    except Exception as e:
        logger.error(f"Error during patched loading: {e}")
        raise
    finally:
        pyimzml.ImzMLParser._get_cv_param = original_get_cv_param # pylint: disable=protected-access
