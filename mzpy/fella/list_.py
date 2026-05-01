"""
Python port of FELLA::listCategories, listMethods, listApprox,
listInternalDatabases
"""

import os
from typing import List


def listCategories() -> List[str]:
    """Return the node categories used in the internal representations."""
    return ["pathway", "module", "enzyme", "reaction", "compound"]


def listMethods() -> List[str]:
    """Return the methods available for the analysis."""
    return ["hypergeom", "diffusion", "pagerank"]


def listApprox() -> List[str]:
    """Return the available approximations for the analysis."""
    return ["simulation", "normality", "gamma", "t"]


def listInternalDatabases(full_names: bool = False) -> List[str]:
    """
    List the directories in the local database path
    (./database relative to the current working directory).
    """
    dir_internal = os.path.join(os.getcwd(), "database")
    if not os.path.isdir(dir_internal):
        print("No local databases have been built yet.")
        return []
    dirs = [
        d for d in os.listdir(dir_internal)
        if os.path.isdir(os.path.join(dir_internal, d))
    ]
    if full_names:
        dirs = [os.path.join(dir_internal, d) for d in dirs]
    return dirs
