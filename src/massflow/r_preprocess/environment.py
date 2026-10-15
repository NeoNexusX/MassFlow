import os
import sys
import subprocess
import shutil
from typing import Optional, Any, Dict
from massflow.tools.logger import get_logger

logger = get_logger("massflow.r_preprocess")

_DEFAULT_R_HOME = None

class RConfig:
    """
    R environment configuration manager.
    Dynamically finds R_HOME and configures system environment variables.
    """
    DEFAULT_R_HOME: Optional[str] = _DEFAULT_R_HOME

    @classmethod
    def detect_system_r_home(cls) -> Optional[str]:
        """
        Get R_HOME path via command line tool.
        """
        if cls.DEFAULT_R_HOME and os.path.exists(cls.DEFAULT_R_HOME):
            return cls.DEFAULT_R_HOME

        if shutil.which("R"):
            output = subprocess.check_output(
                ["R", "RHOME"],
                universal_newlines=True,
                stderr=subprocess.DEVNULL
            )
            return output.strip()

    @staticmethod
    def setup_env(r_home: Optional[str] = None):
        """
        Set environment variables before importing rpy2.
        Priority: Function argument > Environment variable R_HOME > System automatic detection
        """
        # R uses English environment to avoid rpy2 UnicodeDecodeError caused by Chinese errors
        os.environ["LANGUAGE"] = "en"
        os.environ["LC_ALL"] = "C"

        def_vsize = os.environ.get('MASSFLOW_R_MAX_VSIZE', '32Gb')
        def_mem = os.environ.get('MASSFLOW_R_MAX_MEM_SIZE', '32Gb')

        if 'R_MAX_VSIZE' not in os.environ:
            os.environ['R_MAX_VSIZE'] = def_vsize

        if 'R_MAX_MEM_SIZE' not in os.environ:
            os.environ['R_MAX_MEM_SIZE'] = def_mem

        target_home = r_home or os.environ.get("R_HOME") or RConfig.detect_system_r_home()

        if not target_home or not os.path.exists(target_home):
            raise FileNotFoundError(f"R_HOME not found at: {target_home}")

        os.environ['R_HOME'] = target_home

        if sys.platform == 'win32':
            # Windows: Need to add bin/x64 to PATH
            r_bin = os.path.join(target_home, 'bin', 'x64')
            if os.path.exists(r_bin) and r_bin not in os.environ['PATH']:
                os.environ['PATH'] = r_bin + os.pathsep + os.environ['PATH']
            # Force ABI mode to be compatible with different versions
            os.environ.setdefault('RPY2_CFFI_MODE', 'ABI')

        logger.info(f"R Environment configured successfully. R_HOME={target_home}")

class REnvironment:
    """
    Responsible for R environment initialization and lazy loading of rpy2.
    """
    _instance = None
    pkgs: Dict[str, Any] = {}

    def __new__(cls):
        if cls._instance is None:
            instance = super().__new__(cls)
            instance._initialize()
            # Cache only a fully initialized environment so a failed start can be retried.
            cls._instance = instance
        return cls._instance

    def _initialize(self):

        essential_pkgs = ['base', 'utils', 'methods', 'S4Vectors','Cardinal']
        # 1. Configure environment
        try:
            RConfig.setup_env()
            # pylint: disable=import-outside-toplevel
            import rpy2.robjects as robjects
            from rpy2.robjects.packages import importr
            from rpy2.robjects.vectors import FloatVector, IntVector, ListVector
            # pylint: enable=import-outside-toplevel

            # pylint: disable=attribute-defined-outside-init
            self.robjects = robjects
            self.float_vector = FloatVector
            self.int_vector = IntVector
            self.list_vector = ListVector
            self._importr = importr

        # 2. Load base packages
            for pkg in essential_pkgs:
                self.pkgs[pkg] = importr(pkg)
            ver = self.pkgs['utils'].packageVersion('Cardinal')
            logger.info(f"R loaded successfully. Cardinal version: {ver}")

        except ImportError as e:
            logger.error("Required Python package 'rpy2' is not installed.")
            raise ImportError(
                "The Cardinal backend requires rpy2; install it with `pip install massflow[r]`."
            ) from e

        except Exception as e:
            logger.error(f"Failed to initialize R or load the R packages {essential_pkgs}: {e}")
            raise RuntimeError(
                f"Failed to initialize R or load the R packages {essential_pkgs}. "
                "Check R_HOME and that these packages are installed in R."
            ) from e

    @property
    def cardinal(self):
        """Get Cardinal package, reload if namespace is invalid."""
        return self.pkgs.get('Cardinal')
