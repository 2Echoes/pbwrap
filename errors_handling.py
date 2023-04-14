"""This subpackage defines Exceptions classes raised during analysis pipelines."""

class SegmentationError(Exception) :
    pass

class CellnumberError(SegmentationError) :
    pass

class SegmentationProcessError(Exception) :
    pass

class DetectionError(Exception) :
    pass

class PlotError(Exception) :
    pass

class CellanalysisError(Exception) :
    pass
