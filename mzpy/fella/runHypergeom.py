"""
Python port of FELLA::runHypergeom
"""

from typing import Optional

import numpy as np
from scipy import stats

from .AllClasses import FELLA_USER, FELLA_DATA
from .is_ import is_FELLA_USER, is_FELLA_DATA
from .get_ import getStatus, getInput, getBackground, getMatrix, getCom


def runHypergeom(
    object: Optional[FELLA_USER] = None,
    data: Optional[FELLA_DATA] = None,
    p_adjust: str = "fdr",
) -> FELLA_USER:
    """
    Perform over-representation analysis via the hypergeometric test.

    Parameters
    ----------
    object : FELLA_USER
        User object with mapped metabolites.
    data : FELLA_DATA
        Loaded database.
    p_adjust : str
        Method for p-value adjustment (default "fdr").

    Returns
    -------
    FELLA_USER
        Object updated with hypergeometric test results.
    """
    print("Running hypergeom...")

    if not is_FELLA_USER(object):
        print("'object' is not a FELLA_USER object. Returning NULL...")
        return object
    if not is_FELLA_DATA(data):
        print("'data' is not a FELLA_DATA object. Returning NULL...")
        return object
    if getStatus(data) != "loaded":
        print("'data' points to an empty FELLA_DATA object! Returning original object...")
        return object

    print("Starting hypergeometric p-values calculation...")

    hypergeom_matrix = getMatrix(data, "hypergeom")
    if hypergeom_matrix is None or hypergeom_matrix.shape[0] <= 1:
        print("Hypergeometric test failed because its matrix is not loaded.")
        return object

    metabolites_input = getInput(object)
    row_names = data.hypergeom.row_names
    col_names = data.hypergeom.col_names

    metabolites_input_intersect = [c for c in metabolites_input if c in row_names]
    metabolites_background_intersect = [c for c in getBackground(object) if c in row_names]

    if len(metabolites_input) == 0:
        print("Hypergeometric test failed because there are no compounds in the input.")
        return object
    elif len(metabolites_input_intersect) == 0:
        print("None of the compounds were in the hypergeometric test background.")
        return object
    elif len(metabolites_input_intersect) < len(metabolites_input):
        print(
            f"Some compounds have been excluded as they are not in the "
            f"hypergeometric test background. Amount decreased from "
            f"{len(metabolites_input)} to {len(metabolites_input_intersect)}"
        )
        metabolites_input = metabolites_input_intersect

    # If custom background, subset matrix rows
    if len(metabolites_background_intersect) > 0:
        row_idx = [row_names.index(c) for c in metabolites_background_intersect if c in row_names]
        hypergeom_matrix = hypergeom_matrix[row_idx, :]
        row_names = [row_names[i] for i in row_idx]

    row_comp = [i for i, c in enumerate(row_names) if c in metabolites_input]

    n_paths = hypergeom_matrix.shape[1]
    pvalues_path = np.zeros(n_paths)
    pathhits = np.zeros(n_paths)
    pathbackground = np.zeros(n_paths)

    for path in range(n_paths):
        col = hypergeom_matrix[:, path].toarray().flatten().astype(bool)
        sample_success = int(col[row_comp].sum())
        pathhits[path] = sample_success
        total_success = int(col.sum())
        pathbackground[path] = total_success
        total_failure = len(col) - total_success
        sample_size = len(row_comp)

        # phyper(sample_success - 1, total_success, total_failure, sample_size, lower.tail=False)
        pvalues_path[path] = stats.hypergeom.sf(
            sample_success - 1, len(col), total_success, sample_size
        )

    if p_adjust == "fdr":
        from statsmodels.stats.multitest import multipletests
        _, adjusted_pvalues, _, _ = multipletests(pvalues_path, method="fdr_bh")
    else:
        adjusted_pvalues = pvalues_path

    # Store results as pandas Series with names
    import pandas as pd
    object.hypergeom.pvalues = pd.Series(adjusted_pvalues, index=col_names)
    object.hypergeom.nbackground = len(row_names)
    object.hypergeom.ninput = len(metabolites_input_intersect)
    object.hypergeom.pathhits = pd.Series(pathhits, index=col_names)
    object.hypergeom.pathbackground = pd.Series(pathbackground, index=col_names)
    object.hypergeom.valid = True

    print("Done.")
    return object
