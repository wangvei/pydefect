# -*- coding: utf-8 -*-

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"


class StructureError(Exception):
    """Raised when the Structure is inadequate."""
    pass


class CellSizeError(Exception):
    """Raised when the cell size is inadequate."""
    pass


class InvalidFileError(Exception):
    """Raised when the given file is invalid."""
    pass


class NoConvergenceError(Exception):
    """Raised when first-principles calculations are not converged."""
    pass
