"""This subpackage defines Exceptions classes raised during analysis pipelines."""

################# Segmentation Error #####################
class SegmentationError(Exception) :
    """Exception class raised during segmentation."""
    pass

class CellnumberError(SegmentationError) :
    """'SegmentationError' Subclass. Exception raised when segmentation resulted in incorrect cell number."""
    pass

class SegmentationProcessError(SegmentationError) :
    """'SegmentationError' Subclass."""
    pass

class PbodySegmentationError(SegmentationError):
    """'SegmentationError' Subclass."""
    pass
##########################################################
################## Detection Error #######################

class DetectionError(Exception) :
    """Exception class raised during spot detection."""
    pass

class DetectionTimeOutError(DetectionError):
    """'DetectionError' Subclass."""
    pass

class NoSpotError(DetectionError):
    """'DetectionError' Subclass."""
    pass

##########################################################

class PlotError(Exception) :
    """Exception class raised during plots making."""
    pass

class CellanalysisError(Exception) :
    """Exception class raised during cell analysis."""
    pass

class PreprocessingError(Exception) :
    """Exception class raised during preprocessing image."""
    pass