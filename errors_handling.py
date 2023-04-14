"""This subpackage defines Exceptions classes raised during analysis pipelines."""

class SegmentationError(Exception) :
    """Exception class raised during segmentation."""
    pass

class CellnumberError(SegmentationError) :
    """'SegmentationError' Subclass. Exception raised when segmentation resulted in incorrect cell number."""
    pass

class SegmentationProcessError(SegmentationError) :
    """'SegmentationError' Subclass."""
    pass

class DetectionError(Exception) :
    """Exception class raised during spot detection."""
    pass

class DetectionTimeOutError(DetectionError):
    pass

class PlotError(Exception) :
    """Exception class raised during plots making."""
    pass

class CellanalysisError(Exception) :
    """Exception class raised during cell analysis."""
    pass

class PreprocessingError(Exception) :
    """Exception class raised during preprocessing image."""
    pass