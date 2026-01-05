from massflow.module.mass_spectrum_set import MassSpectrumSet
from massflow.module.ms_data_manager_imzml import MSDataManagerImzML
from massflow.r_preprocess.environment import REnvironment
from massflow.logger import get_logger

logger = get_logger("r_preprocess_adapter")

class CardinalAdapter:
    """
    Preprocessing adapter, calling the R language Cardinal package for data preprocessing.
    """
    @staticmethod
    def align_massdata(dm_data: MSDataManagerImzML,
                       reference=None,
                       tolerance=None,
                       units: str = "ppm",
                       binfun: str = "min",
                       binratio: float = 2.0,
                       temp_dir: str = "./temp_align_data"
                       ) -> MSDataManagerImzML:
        """
        Call Cardinal::peakAlign for peak alignment.
        """
        r_env = REnvironment()

        ms_cardinal = MassSpectrumSet()
        dm_cardinal = MSDataManagerImzML(ms_cardinal, temp_dir=temp_dir)
        dm_cardinal.copy_meta(dm_data)
        dm_cardinal.ms.meta.continuous = True

        r_reference = r_env.FloatVector(reference) if reference is not None else r_env.robjects.NULL
        r_tol = r_env.FloatVector([tolerance]) if tolerance is not None else r_env.robjects.NA_Real

        logger.info(f"Starting peak alignment using Cardinal::peakAlign")

        imzml_filepath = dm_data.filepath
        r_massdata = r_env.cardinal.readImzML(imzml_filepath)

        aligned_massdata = r_env.cardinal.peakAlign(
            r_massdata,
            ref=r_reference,
            tolerance=r_tol,
            units=units,
            binfun=binfun,
            binratio=binratio
        )

        # 5. Write results
        aligned_filepath = dm_cardinal.filepath
        r_env.cardinal.writeMSIData(aligned_massdata, file=aligned_filepath, bundle=False)


        dm_cardinal.load_full_data_from_file()

        logger.info(f"Peak alignment completed and data saved to {aligned_filepath}")
        return dm_cardinal
    
    @staticmethod
    def peak_pick(dm_data: MSDataManagerImzML,
                  method: str = "mad",
                  SNR: float = 3.0,
                  temp_dir: str = "./temp_pick_data"
                  ) -> MSDataManagerImzML:
        """
        Call Cardinal::peakPick for peak picking.
        """
        r_env = REnvironment()
        
        ms_cardinal = MassSpectrumSet()
        dm_cardinal = MSDataManagerImzML(ms_cardinal, temp_dir=temp_dir)
        dm_cardinal.copy_meta(dm_data)
        dm_cardinal.ms.meta.processed = True
        

        r_snr = r_env.FloatVector([SNR])

        logger.info(f"Starting peak picking using Cardinal::peakPick (method={method}, SNR={SNR})")

        imzml_filepath = dm_data.filepath
        r_massdata = r_env.cardinal.readImzML(imzml_filepath)

        picked_massdata = r_env.cardinal.peakPick(
            r_massdata,
            method=method,
            SNR=r_snr
        )

        picked_massdata_realize = r_env.cardinal.process(picked_massdata)

        picked_filepath = dm_cardinal.filepath
        r_env.cardinal.writeMSIData(picked_massdata_realize, file=picked_filepath, bundle=False)
        
        dm_cardinal.load_full_data_from_file()

        logger.info(f"Peak picking completed and data saved to {picked_filepath}")
        return dm_cardinal