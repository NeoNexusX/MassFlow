from massflow.module import MassSpectrumSet
from massflow.data_manager import MSDataManagerImzML
from massflow.r_preprocess.environment import REnvironment
from massflow.tools.imzml_monkey_patch import risky_imzml_loader
from massflow.tools.logger import get_logger

logger = get_logger("massflow.r_preprocess")

class CardinalAdapter:
    """
    Preprocessing adapter, calling the R language Cardinal package for data preprocessing.
    """
    @staticmethod
    def peak_align(data_manager: MSDataManagerImzML,
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
        # Initialize shared R runtime and Cardinal bindings.
        r_env = REnvironment()

        # Build output data manager and copy metadata from source dataset.
        ms_cardinal = MassSpectrumSet()
        dm_cardinal = MSDataManagerImzML(ms_cardinal, temp_dir=temp_dir)
        dm_cardinal.copy_meta(data_manager)
        dm_cardinal.ms.meta.continuous = True

        # Convert optional Python parameters to R-compatible values.
        # - reference=None -> NULL (let Cardinal auto-select reference)
        # - tolerance=None -> NA_real_ (use package default tolerance)
        r_reference = r_env.float_vector(reference) if reference is not None else r_env.robjects.NULL
        r_tol = r_env.float_vector([tolerance]) if tolerance is not None else r_env.robjects.NA_Real

        logger.info(f"Starting peak alignment using Cardinal::peakAlign(units={units}, binfun={binfun}, binratio={binratio})")

        # Read source imzML as Cardinal MSImagingExperiment.
        imzml_filepath = data_manager.imzml_filepath
        r_massdata = r_env.cardinal.readImzML(imzml_filepath) # type: ignore

        # Run Cardinal peak alignment in R.
        aligned_massdata = r_env.cardinal.peakAlign( # type: ignore
            r_massdata,
            ref=r_reference,
            tolerance=r_tol,
            units=units,
            binfun=binfun,
            binratio=binratio
        )

        # Persist aligned result to output imzML and reload Python-side headers.
        aligned_filepath = dm_cardinal.imzml_filepath
        r_env.cardinal.writeMSIData(aligned_massdata, file=aligned_filepath, bundle=False) # type: ignore


        with risky_imzml_loader():
            dm_cardinal.load_head_data()

        logger.info(f"Peak alignment completed and data saved to {aligned_filepath}")
        return dm_cardinal

    @staticmethod
    def peak_pick(data_manager: MSDataManagerImzML,
                  width: int = 5,
                  method: str = "quantile", # "diff", "sd", "mad", "quantile", "filter", "cwt"
                  snr: float = 2.0,
                  return_type: str = "height", # "height", "area"
                  prominence: float | None = None,
                  relheight: float | None = None,
                  nbins: int = 1,
                  overlap: float = 0.5,
                  temp_dir: str = "./temp_pick_data"
                  ) -> MSDataManagerImzML:
        """
        Call Cardinal::peakPick for peak picking.
        """
        if relheight is not None:
            logger.warning("The 'relheight' parameter is not supported in the Cardinal backend and will be ignored.")

        # Initialize shared R runtime and Cardinal bindings.
        r_env = REnvironment()

        # Build output data manager and copy metadata from source dataset.
        ms_cardinal = MassSpectrumSet()
        dm_cardinal = MSDataManagerImzML(ms_cardinal, temp_dir=temp_dir)
        dm_cardinal.copy_meta(data_manager)

        logger.info(f"Starting peak picking using Cardinal::peakPick (method={method}, SNR={snr}, type={return_type})")

        # Read source imzML as Cardinal MSImagingExperiment.
        imzml_filepath = data_manager.imzml_filepath
        r_massdata = r_env.cardinal.readImzML(imzml_filepath) # type: ignore
        r_as = r_env.robjects.r["as"]

        r_massdata_arrays = r_as(r_massdata, "MSImagingArrays") # type: ignore

        # Run Cardinal peak picking in R.
        picked_massdata = r_env.cardinal.peakPick( # type: ignore
            r_massdata_arrays,
            width=width,
            method=method,
            SNR=snr,
            type=return_type,
            prominence=prominence if prominence is not None else r_env.robjects.NULL,
            nbins=nbins,
            overlap=overlap
        )

        # Materialize the queued peak-pick on the arrays representation before export.
        picked_massdata_realize = r_env.cardinal.process(picked_massdata) # type: ignore

        # Persist picked result to output imzML and reload Python-side headers.
        picked_filepath = dm_cardinal.imzml_filepath
        r_env.cardinal.writeMSIData(picked_massdata_realize, file=picked_filepath, bundle=False) # type: ignore

        with risky_imzml_loader():
            dm_cardinal.load_head_data()

        logger.info(f"Peak picking completed and data saved to {picked_filepath}")
        return dm_cardinal
