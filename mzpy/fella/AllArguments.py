"""
Python port of FELLA:::checkArguments

Performs argument validation mimicking the original R function.
Each parameter is only checked if it is passed (not None).
"""

from typing import Optional, List, Any

from .AllClasses import FELLA_DATA, FELLA_USER
from .list_ import listMethods, listApprox


def checkArguments(
    databaseDir: Optional[str] = None,
    internalDir: Optional[bool] = None,
    method: Optional[str] = None,
    methods: Optional[List[str]] = None,
    approx: Optional[str] = None,
    loadMatrix: Optional[List[str]] = None,
    threshold: Optional[float] = None,
    plimit: Optional[int] = None,
    nlimit: Optional[int] = None,
    niter: Optional[int] = None,
    layout: Optional[bool] = None,
    thresholdConnectedComponent: Optional[float] = None,
    dampingFactor: Optional[float] = None,
    t_df: Optional[float] = None,
    object: Optional[Any] = None,
    data: Optional[Any] = None,
    **kwargs,
) -> dict:
    """
    Check arguments and return {'valid': True/False, 'ans': ...}.
    Mimics the original R behaviour.
    """
    # databaseDir
    if databaseDir is not None:
        if not isinstance(databaseDir, str) or len(databaseDir) == 0:
            print("'databaseDir' must be a non-empty string.")
            return {"valid": False, "ans": None}

    # internalDir
    if internalDir is not None:
        if not isinstance(internalDir, bool):
            print("'internalDir' must be a boolean.")
            return {"valid": False, "ans": None}

    # method
    if method is not None:
        if not isinstance(method, str) or method not in listMethods():
            print(
                "'method' must be one of: 'hypergeom', 'diffusion', 'pagerank'."
            )
            return {"valid": False, "ans": object}

    # methods
    if methods is not None:
        if not isinstance(methods, (list, tuple)) or not all(
            m in listMethods() for m in methods
        ):
            print(
                "'methods' must be a list containing 'hypergeom', "
                "'diffusion' and/or 'pagerank'."
            )
            return {"valid": False, "ans": object}

    # approx
    if approx is not None:
        if approx not in listApprox():
            print(
                "'approx' must be one of: 'simulation', 'normality', "
                "'gamma', 't'."
            )
            return {"valid": False, "ans": None}

    # loadMatrix
    if loadMatrix is not None:
        if not isinstance(loadMatrix, (list, tuple)):
            print("'loadMatrix' must be a list/tuple or None.")
            return {"valid": False, "ans": None}
        valid_lm = {"diffusion", "pagerank"}
        if not set(loadMatrix).issubset(valid_lm):
            print(
                "'loadMatrix' can only contain 'diffusion' and/or 'pagerank'."
            )
            return {"valid": False, "ans": None}

    # threshold
    if threshold is not None:
        if not isinstance(threshold, (int, float)) or not (0 < threshold <= 1):
            print("'threshold' must be numeric between 0 and 1.")
            return {"valid": False, "ans": None}

    # plimit
    if plimit is not None:
        if not isinstance(plimit, int) or not (1 <= plimit <= 50):
            print("'plimit' must be an integer between 1 and 50.")
            return {"valid": False, "ans": None}

    # nlimit
    if nlimit is not None:
        if not isinstance(nlimit, int) or not (1 <= nlimit <= 1000):
            print("'nlimit' must be an integer between 1 and 1000.")
            return {"valid": False, "ans": None}

    # niter
    if niter is not None:
        if not isinstance(niter, int) or not (100 <= niter <= 100000):
            print("'niter' must be an integer between 1e2 and 1e5.")
            return {"valid": False, "ans": None}

    # layout
    if layout is not None:
        if not isinstance(layout, bool):
            print("'layout' must be a boolean.")
            return {"valid": False, "ans": None}

    # thresholdConnectedComponent
    if thresholdConnectedComponent is not None:
        if (
            not isinstance(thresholdConnectedComponent, (int, float))
            or not (0 <= thresholdConnectedComponent <= 1)
        ):
            print(
                "'thresholdConnectedComponent' must be numeric between 0 and 1."
            )
            return {"valid": False, "ans": None}

    # dampingFactor
    if dampingFactor is not None:
        if (
            not isinstance(dampingFactor, (int, float))
            or not (0 < dampingFactor < 1)
        ):
            print(
                "'dampingFactor' must be numeric between 0 and 1 (exclusive)."
            )
            return {"valid": False, "ans": None}

    # t.df
    if t_df is not None:
        if not isinstance(t_df, (int, float)) or t_df <= 0:
            print("'t.df' must be a positive numeric value.")
            return {"valid": False, "ans": None}

    # object
    if object is not None:
        if not isinstance(object, FELLA_USER):
            print("'object' is not a FELLA_USER object.")
            return {"valid": False, "ans": None}

    # data
    if data is not None:
        if not isinstance(data, FELLA_DATA):
            print("'data' is not a FELLA_DATA object.")
            return {"valid": False, "ans": None}

    return {"valid": True, "ans": None}
