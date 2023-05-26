class SetError(Exception) :
    """Error raise solutions are outside of the given Set."""
    pass
class SolutionNotRealError(SetError):
    """Error raised when solution is not real"""
    pass
