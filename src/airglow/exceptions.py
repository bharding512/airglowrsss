"""Custom exceptions for the FPI data processing pipeline."""

class ProcessingError(Exception):
    """Base exception for all processing-related errors"""
    pass

class LaserProcessingError(ProcessingError):
    """Base exception for all laser-related processing errors"""
    pass

class BadLaserError(LaserProcessingError):
    """Raised when laser centerfinding fails"""
    pass

class InstrumentProcessingError(ProcessingError):
    """Base exception for general instrument processing errors"""
    pass

class BadDopplerReferenceError(ProcessingError):
    """Raised when a Doppler reference cannot be established"""
    pass

class NoSkyImagesError(ProcessingError):
    """Raised when no sky images are found for a given sky line tag"""
    pass