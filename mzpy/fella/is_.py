"""
Python port of FELLA:::is.FELLA.DATA and is.FELLA.USER
"""

from .AllClasses import FELLA_DATA, FELLA_USER


def is_FELLA_DATA(x) -> bool:
    """Check if x is a FELLA_DATA object."""
    return isinstance(x, FELLA_DATA)


def is_FELLA_USER(x) -> bool:
    """Check if x is a FELLA_USER object."""
    return isinstance(x, FELLA_USER)
