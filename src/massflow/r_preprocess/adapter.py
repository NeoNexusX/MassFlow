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
    _r_env: REnvironment | None = None

    @classmethod
    def r_env(cls) -> REnvironment:
        if cls._r_env is None:
            cls._r_env = REnvironment()
        return cls._r_env

    @classmethod
    def peak_align(
        cls,
        data_manager: MSDataManagerImzML,
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
        r_env = cls.r_env()
        r_reference = r_env.float_vector(reference) if reference is not None else r_env.robjects.NULL
        r_tolerance = r_env.float_vector([tolerance]) if tolerance is not None else r_env.robjects.NA_Real

        logger.info(f"Starting peak alignment using Cardinal::peakAlign(units={units}, binfun={binfun}, binratio={binratio})")
        return cls.cardinal_runner(
            data_manager=data_manager,
            func="peakAlign",
            temp_dir=temp_dir,
            ref=r_reference,
            tolerance=r_tolerance,
            units=units,
            binfun=binfun,
            binratio=binratio,
        )

    @classmethod
    def peak_pick(
        cls,
        data_manager: MSDataManagerImzML,
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

        r_env = cls.r_env()
        r_prominence = prominence if prominence is not None else r_env.robjects.NULL

        logger.info(f"Starting peak picking using Cardinal::peakPick (method={method}, SNR={snr}, type={return_type})")
        return cls.cardinal_runner(
            data_manager=data_manager,
            func="peakPick",
            temp_dir=temp_dir,
            realize=True,
            width=width,
            method=method,
            SNR=snr,
            type=return_type,
            prominence=r_prominence,
            nbins=nbins,
            overlap=overlap,
        )

    @classmethod
    def baseline_correction(
        cls,
        data_manager: MSDataManagerImzML,
        method: str = "locmin",
        smooth: str = "none",
        span: float = 0.1,
        upper: bool = False,
        width: int = 5,
        decreasing: bool = True,
        temp_dir: str = "./temp_baseline_data",
    ) -> MSDataManagerImzML:
        """
        Call Cardinal::reduceBaseline for baseline correction.
        """
        cardinal_method = method.strip().lower().removesuffix("_numba")
        logger.info(f"Starting baseline reduction using Cardinal::reduceBaseline(method={cardinal_method})")

        kwargs = {"method": cardinal_method}
        if cardinal_method == "locmin":
            kwargs.update({ # type: ignore
                "smooth": smooth,
                "span": span,
                "upper": upper,
            })
        elif cardinal_method == "snip":
            kwargs["width"] = int(width) # type: ignore
            kwargs["decreasing"] = decreasing # type: ignore
        elif cardinal_method in {"median", "med"}:
            kwargs["method"] = "median"
            kwargs["width"] = int(width) # type: ignore
        elif cardinal_method == "hull":
            kwargs["upper"] = upper # type: ignore

        return cls.cardinal_runner(
            data_manager=data_manager,
            func="reduceBaseline",
            temp_dir=temp_dir,
            realize=True,
            **kwargs,
        )

    @classmethod
    def noise_reduction(
        cls,
        data_manager: MSDataManagerImzML,
        method: str = "ma",
        window: int = 5,
        sd: float | None = None,
        polyorder: int = 3,
        deriv: int = 0,
        delta: float = 1.0,
        temp_dir: str = "./temp_noise_data",
    ) -> MSDataManagerImzML:
        """
        Call Cardinal::smooth for spectral smoothing/noise reduction.
        """
        method_norm = method.strip().lower()
        cardinal_method = "sgolay" if method_norm in {"savgol", "savgol_numba"} else method_norm
        cardinal_method = cardinal_method.removesuffix("_numba")

        logger.info(f"Starting noise reduction using Cardinal::smooth(method={cardinal_method})")

        kwargs = {
            "method": cardinal_method,
            "width": int(window),
        }
        if cardinal_method == "gaussian" and sd is not None:
            kwargs["sd"] = sd
        elif cardinal_method == "sgolay":
            kwargs.update({
                "order": int(polyorder),
                "deriv": int(deriv),
                "delta": float(delta),
            })

        return cls.cardinal_runner(
            data_manager=data_manager,
            func="smooth",
            temp_dir=temp_dir,
            realize=True,
            **kwargs,
        )

    @classmethod
    def normalization(
        cls,
        data_manager: MSDataManagerImzML,
        method: str = "tic",
        scale: float | None = None,
        ref: float | None = None,
        temp_dir: str = "./temp_normalization_data",
    ) -> MSDataManagerImzML:
        """
        Call Cardinal::normalize for intensity normalization.
        """
        r_env = cls.r_env()
        method_norm = method.strip().lower().removesuffix("_numba")
        cardinal_method = "reference" if method_norm in {"ref", "reference"} else method_norm

        logger.info(f"Starting normalization using Cardinal::normalize(method={cardinal_method})")

        kwargs = {"method": cardinal_method}
        if scale is not None:
            kwargs["scale"] = float(scale) # type: ignore
        if cardinal_method == "reference":
            if ref is None:
                raise ValueError("Cardinal reference normalization requires ref")
            kwargs["ref"] = r_env.float_vector([ref]) # type: ignore

        return cls.cardinal_runner(
            data_manager=data_manager,
            func="normalize",
            temp_dir=temp_dir,
            realize=True,
            **kwargs,
        )

    @classmethod
    def normalize(cls, **kwargs) -> MSDataManagerImzML:
        """
        Alias for normalization to match the Cardinal function name.
        """
        return cls.normalization(**kwargs)

    @classmethod
    def cardinal_runner(
        cls,
        data_manager: MSDataManagerImzML,
        func: str,
        temp_dir: str = "./temp_cardinal",
        *,
        realize: bool = False,
        **kwargs,
    ) -> MSDataManagerImzML:
        """
        Run a Cardinal R function for preprocessing.
        """
        # Initialize shared R runtime and Cardinal bindings.
        r_env = cls.r_env()
        cardinal_func = getattr(r_env.cardinal, func, None)
        if cardinal_func is None:
            try:
                cardinal_func = r_env.robjects.r[func]
            except LookupError:
                cardinal_func = None
        if cardinal_func is None:
            raise AttributeError(f"Cardinal function '{func}' is not available.")

        # Build output data manager and copy metadata from source dataset.
        dm_cardinal = MSDataManagerImzML(MassSpectrumSet(), temp_dir=temp_dir)
        dm_cardinal.copy_meta(data_manager)

        # Read source imzML as Cardinal MSImagingExperiment.
        input_filepath = data_manager.imzml_filepath
        r_massdata = r_env.cardinal.readImzML(input_filepath) # type: ignore
        r_input = r_massdata

        if func == "peakPick":
            r_as = r_env.robjects.r["as"]
            r_input = r_as(r_massdata, "MSImagingArrays") # type: ignore

        result_massdata = cardinal_func(r_input, **kwargs) # type: ignore

        if realize:
            result_massdata = r_env.cardinal.process(result_massdata) # type: ignore

        # Persist processed result to output imzML and reload Python-side headers.
        output_filepath = dm_cardinal.imzml_filepath
        r_env.cardinal.writeMSIData(result_massdata, file=output_filepath, bundle=False) # type: ignore

        with risky_imzml_loader():
            dm_cardinal.load_head_data()

        logger.info(f"Cardinal::{func} completed and data saved to {output_filepath}")
        return dm_cardinal
