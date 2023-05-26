################# Segmentation Error #####################
class SegmentationError(Exception) :
    """Exception class raised during segmentation."""
    pass

class CellSegmentationError(SegmentationError) :
    """SegmentationError Subclass. Genereal Error raised when an expection occurs during CellSegmentation """

class CellnumberError(CellSegmentationError) :
    """SegmentationError' Subclass. Exception raised when segmentation resulted in incorrect cell number."""
    pass

class PbodySegmentationError(SegmentationError):
    """SegmentationError Subclass. Genereal Error raised when an expection occurs during P-bodies segmentation"""
    pass