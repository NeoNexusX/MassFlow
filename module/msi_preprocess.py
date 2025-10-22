"""
MSI Preprocessing Module

Provides abstract base class and concrete implementations for MSI data preprocessing.
Includes peak picking, TIC normalization, peak alignment, baseline correction, and noise reduction.
"""
import numpy as np
from .msi_module import MSI

class MSIPreprocessor():
    """
    Abstract base class for MSI data preprocessing.
    
    This class provides a framework for implementing various preprocessing techniques
    for Mass Spectrometry Imaging (MSI) data, including peak picking, normalization,
    alignment, baseline correction, and noise reduction.
    
    Supports both traditional input (mz, msroi arrays) and MSI object input for
    better integration with the MSI framework.
    
    Attributes:
        msi_object (MSI, optional): MSI object containing data and metadata
        preprocessing_params (dict): Parameters for preprocessing operations
        processed_data (np.ndarray, optional): Processed MSI data
    """

    def __init__(self, msi: MSI):
        """
        Initialize MSI preprocessor.
        
        Args:
            msi_object (MSI, optional): MSI object containing data and metadata
            preprocessing_params (dict, optional): Parameters for preprocessing operations
            
        Raises:
            ValueError: If neither msi_object nor (mz, msroi) are provided
        """
        self.processed_data = msi.copy()
        self.original_data = msi

    def peak_pick(self,method: str):
        """
        Perform peak picking on MSI data.
        
        Abstract method that must be implemented by concrete classes to identify
        significant peaks in the mass spectra.
        
        Args:
            Fill this
            
        Returns:
            processed_data (MSI): Processed MSI object with peak-picked data
            
        Raises:
            NotImplementedError: If not implemented by subclass
        """

    def tic_normalization(self,method: str = 'total_ion_current') -> MSI:
        """
        Perform Total Ion Current (TIC) normalization.
        
        Abstract method for normalizing MSI data using various TIC-based approaches
        to account for variations in total ion intensity across pixels.
        
        Args:
            method (str): Normalization method ('total_ion_current', 'median', 'max')
            
        Returns:
            processed_data (MSI): Processed MSI object with peak-picked data
            
        Raises:
            NotImplementedError: If not implemented by subclass
        """

    def peak_alignment(self)-> MSI:
        """
        Perform peak alignment across spectra.
        
        Abstract method for aligning peaks across different spectra to correct
        for mass calibration drift and improve data consistency.
        
        Args:
            fill this
            
        Returns:
            processed_data (MSI): Processed MSI object with peak_aligned data
            
        Raises:
            NotImplementedError: If not implemented by subclass
        """

    def baseline_correction(self) -> MSI:
        """
        Perform baseline correction on MSI data.
        
        Abstract method for removing baseline drift and background signals
        from mass spectra to improve peak detection and quantification.
        
        Args:
            fill this
            
        Returns:
            processed_data (MSI): Processed MSI object with baseline_correction
            
        Raises:
            NotImplementedError: If not implemented by subclass
        """

    def noise_reduction(self) -> MSI:
        """
        Perform noise reduction on MSI data.
        
        Abstract method for reducing noise in mass spectra while preserving
        important spectral features and peak information.
        
        Args:
            fill this
        Returns:
            processed_data (MSI): Processed MSI object with noise_reduction
            
        Raises:
            NotImplementedError: If not implemented by subclass
        """

    def preprocess_pipeline(self) -> np.ndarray:
        """
        Execute a complete preprocessing pipeline.
        
        Applies multiple preprocessing steps in sequence to the MSI data.
        
        Args:
            steps (List[str]): List of preprocessing steps to apply
            **kwargs: Additional parameters for preprocessing methods
            
        Returns:
            np.ndarray: Fully preprocessed MSI data
            
        Raises:
            ValueError: If invalid preprocessing step is specified
        """
        return self.processed_data
