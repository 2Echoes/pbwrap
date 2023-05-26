class DetectionError(Exception) :
    """Exception class raised during spot detection."""
    pass

class DetectionTimeOutError(DetectionError):
    """'DetectionError' Subclass raised when detection takes too much time."""
    pass

class NoSpotError(DetectionError):
    """'DetectionError' Subclass raised when no spots or too few spots have been detected."""
    pass

class TooManySpotsError(DetectionError):
    """'DetectionError' Subclass raised when too much spots have been detected."""
    pass
