"""
Python port of FELLA::defineCompounds
"""

from typing import List, Optional

from .AllClasses import FELLA_USER, U_userinput
from .is_ import is_FELLA_DATA
from .get_ import getCom, getStatus


def defineCompounds(
    compounds: Optional[List[str]] = None,
    compoundsBackground: Optional[List[str]] = None,
    data = None,
) -> FELLA_USER:
    """
    Create a FELLA_USER object from a list of compounds and a FELLA_DATA object.

    Parameters
    ----------
    compounds : list of str or None
        KEGG compound IDs to enrich.
    compoundsBackground : list of str or None
        Custom background compound IDs.
    data : FELLA_DATA
        Loaded FELLA_DATA object.

    Returns
    -------
    FELLA_USER
        Object with mapped metabolites ready for enrichment.
    """
    if not is_FELLA_DATA(data):
        raise ValueError("'data' is not a FELLA_DATA object")
    if getStatus(data) != "loaded":
        raise ValueError("'data' points to an empty FELLA_DATA object")

    obj = FELLA_USER()

    # background
    if compoundsBackground is None or len(compoundsBackground) == 0:
        print("No background compounds specified. Default background will be used.")
        compoundsBackground = getCom(data, "compound")
    else:
        compoundsBackground = list(set(str(c) for c in compoundsBackground))
        compoundsBackground = [
            c for c in compoundsBackground if c in getCom(data, "compound")
        ]
        if len(compoundsBackground) < 10:
            print(
                "Less than ten specified background compounds appear in KEGG. "
                "Default background will be used instead."
            )
            compoundsBackground = getCom(data, "compound")
        else:
            obj.userinput.metabolitesbackground = compoundsBackground

    # input compounds
    if compounds is None or len(compounds) == 0:
        raise ValueError("Argument 'compounds' cannot be empty.")
    compounds = list(set(str(c) for c in compounds))
    compounds_old = compounds.copy()

    compoundsInBackground = [c for c in compounds if c in compoundsBackground]
    if len(compoundsInBackground) < len(compounds):
        print(
            "Some compounds were introduced as affected but do not belong to "
            "the background. These compounds will be excluded from the analysis. "
            "Use 'getExcluded' to see them."
        )
        compounds = compoundsInBackground

    if len(compounds) == 0:
        raise ValueError(
            "None of the specified compounds appear in the available KEGG data."
        )
    else:
        obj.userinput.metabolites = compounds

    obj.userinput.excluded = [c for c in compounds_old if c not in compounds]
    return obj
