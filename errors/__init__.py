# -*- coding: utf-8 -*-
# Author: Floric Slimani <floric.slimani@live.fr>


"""
This pbwrap subpackage contains custom Exception class for pbwrap usage.
"""
from .segmentation import CellnumberError,SegmentationError,CellSegmentationError,PbodySegmentationError

from .CurveAnalysis import SetError,SolutionNotRealError

from .detection import NoSpotError,DetectionError,TooManySpotsError,DetectionTimeOutError

from .other import PlotError,CellanalysisError,PreprocessingError