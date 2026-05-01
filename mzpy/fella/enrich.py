"""
Python port of FELLA::enrich

Wrapper to perform enrichment analysis in a single call.
"""

from typing import List, Optional, Any

from .AllClasses import FELLA_DATA, FELLA_USER
from .is_ import is_FELLA_DATA
from .loadKEGGdata import loadKEGGdata
from .defineCompounds import defineCompounds
from .runHypergeom import runHypergeom
from .runDiffusion import runDiffusion
from .runPagerank import runPagerank
from .list_ import listMethods


def enrich(
    compounds: Optional[List[str]] = None,
    compoundsBackground: Optional[List[str]] = None,
    methods: Optional[List[str]] = None,
    loadMatrix: str = "none",
    approx: str = "normality",
    t_df: float = 10.0,
    niter: int = 1000,
    h5_path: Optional[str] = None,
    data: Optional[FELLA_DATA] = None,
    **kwargs,
) -> Any:
    """
    Enrich a list of metabolites.

    Parameters
    ----------
    compounds : list of str or None
        KEGG compound IDs to enrich.
    compoundsBackground : list of str or None
        Custom background compound IDs.
    methods : list of str or None
        Methods to run (default all).
    loadMatrix : str
        Deprecated compatibility argument; use loadKEGGdata directly if needed.
    approx : str
        Parametric approximation for diffusion/pagerank.
    t_df : float
        Degrees of freedom for t approximation.
    niter : int
        Permutations for simulation.
    h5_path : str or None
        Path to the HDF5 database file. Required if ``data`` is not supplied.
    data : FELLA_DATA or None
        Pre-loaded database.
    **kwargs
        Extra arguments passed to runDiffusion/runPagerank.

    Returns
    -------
    FELLA_USER or dict
        If data was supplied, returns the updated FELLA_USER.
        Otherwise returns a dict {'user': FELLA_USER, 'data': FELLA_DATA}.
    """
    if methods is None:
        methods = listMethods()

    return_list = False
    if not is_FELLA_DATA(data):
        if h5_path is None:
            raise ValueError("Either 'data' or 'h5_path' must be provided.")
        print("No data object supplied. Loading it from the 'h5_path' argument...")
        return_list = True
        # loadMatrix argument in original can be "none" or c("diffusion","pagerank")
        lm = None if loadMatrix == "none" else [loadMatrix] if isinstance(loadMatrix, str) else list(loadMatrix)
        data = loadKEGGdata(
            h5_path=h5_path,
            load_matrix=lm,
        )

    obj = defineCompounds(
        compounds=compounds,
        compoundsBackground=compoundsBackground,
        data=data,
    )

    if "hypergeom" in methods:
        obj = runHypergeom(object=obj, data=data)
    if "diffusion" in methods:
        obj = runDiffusion(
            object=obj, data=data, approx=approx, t_df=t_df, niter=niter, **kwargs
        )
    if "pagerank" in methods:
        obj = runPagerank(
            object=obj, data=data, approx=approx, dampingFactor=kwargs.get("dampingFactor", 0.85), t_df=t_df, niter=niter, **kwargs
        )

    if return_list:
        return {"user": obj, "data": data}
    return obj
