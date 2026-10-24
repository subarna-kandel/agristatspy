"""
AgriStatsPy custom exceptions.

All package-specific errors inherit from AgriStatsError so that callers
can catch the entire hierarchy with a single except clause.
"""


class AgriStatsError(Exception):
    """Base exception for all AgriStatsPy errors."""


class DesignError(AgriStatsError):
    """Raised when design generation fails due to invalid parameters."""


class ValidationError(AgriStatsError):
    """Raised when a generated design fails its internal invariant checks."""


class AnalysisError(AgriStatsError):
    """Raised when statistical analysis cannot be performed on the supplied data."""


class ExportError(AgriStatsError):
    """Raised when exporting a design or result fails."""
