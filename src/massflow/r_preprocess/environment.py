import os
import sys
import subprocess
import shutil
from typing import Optional, Any, Dict
from massflow.logger import get_logger

logger = get_logger("R Environment")

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
            try:
                output = subprocess.check_output(
                    ["R", "RHOME"], 
                    universal_newlines=True, 
                    stderr=subprocess.DEVNULL
                )
                return output.strip()
            except subprocess.CalledProcessError:
                pass
        
        logger.error("R_HOME could not be detected automatically. Please set R_HOME environment variable or provide R_HOME path.")
        return None
        
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
    
    # rpy2 object cache
    robjects = None
    FloatVector = None
    IntVector = None
    ListVector = None

    def __new__(cls):
        if cls._instance is None:
            cls._instance = super(REnvironment, cls).__new__(cls)
            cls._instance._initialize()
        return cls._instance

    def _initialize(self):
        try:
            # 1. Configure environment
            RConfig.setup_env()
            
            # 2. Lazy import
            try:
                import rpy2.robjects as robjects
                from rpy2.robjects.packages import importr
                from rpy2.robjects.vectors import FloatVector, IntVector, ListVector
            except ImportError as e:
                logger.error("Required Python package 'rpy2' is not installed.")
                raise ImportError("Please install rpy2 via 'pip install rpy2' before running.") from e

            self.robjects = robjects
            self.FloatVector = FloatVector
            self.IntVector = IntVector
            self.ListVector = ListVector
            self._importr = importr

            # 3. Load base packages
            try:
                self.pkgs['base'] = importr('base')
                self.pkgs['utils'] = importr('utils')
                self.pkgs['methods'] = importr('methods')
                self.pkgs['s4'] = importr('S4Vectors')
            except Exception as e:
                logger.error("Failed to load essential R packages.")
                raise RuntimeError("Could not load base R packages.") from e
            
            # 4. Try to load Cardinal
            try:
                self.pkgs['cardinal'] = importr('Cardinal')
            except Exception as e:
                logger.error("Failed to load 'Cardinal' R package. Ensure it is installed in your R environment.")
                raise RuntimeError("Could not load Cardinal package.") from e
            
            ver = self.pkgs['utils'].packageVersion('Cardinal')
            logger.info(f"R loaded successfully. Cardinal version: {ver}")
        
        except (ImportError, RuntimeError):
            raise

        except Exception as e:
            logger.error(f"Unexpected R Environment initialization failure: {e}")
            raise RuntimeError(f"Unknown error during R setup: {type(e).__name__}") from e
        
    @property
    def cardinal(self):
        """Get Cardinal package, reload if namespace is invalid."""
        try:
            # Try to access a Cardinal function to verify if the package is available
            _ = self.pkgs['cardinal'].readImzML
            return self.pkgs['cardinal']
        except Exception:
            # Reload Cardinal package
            logger.warning("Cardinal package namespace invalid, reloading...")
            self.pkgs['cardinal'] = self._importr('Cardinal')
            return self.pkgs['cardinal']