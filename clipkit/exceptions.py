class ClipKITException(Exception):
    pass


class InvalidInputFileFormat(ClipKITException):
    pass


class StopCodonValidationError(ValueError):
    pass
