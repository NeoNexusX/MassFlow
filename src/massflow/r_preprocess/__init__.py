import os
from typing import Optional
from massflow.r_preprocess.environment import RConfig, REnvironment

def set_default_r_home(r_home: str):
    """
    Set the default R_HOME path.
    """
    RConfig.DEFAULT_R_HOME = r_home

    if not r_home or not os.path.exists(r_home):
        raise FileNotFoundError(f"R_HOME not found at: {r_home}")
    
def init_r_environment(r_home: Optional[str] = None):
    """
    Initialize R environment by setting up necessary environment variables.
    """
    if r_home:
        RConfig.setup_env(r_home=r_home)
    return REnvironment()