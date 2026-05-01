"""
Python port of FELLA::runDiffusion
"""

from typing import Optional

import numpy as np
from scipy import stats
from scipy import linalg as la

from .AllClasses import FELLA_USER, FELLA_DATA
from .is_ import is_FELLA_USER, is_FELLA_DATA
from .get_ import (
    getStatus, getInput, getBackground, getMatrix, getCom,
    getGraph, getSums,
)


def runDiffusion(
    object: Optional[FELLA_USER] = None,
    data: Optional[FELLA_DATA] = None,
    approx: str = "normality",
    t_df: float = 10.0,
    niter: int = 1000,
) -> FELLA_USER:
    """
    Perform diffusion-based enrichment.

    Parameters
    ----------
    object : FELLA_USER
        User object with mapped metabolites.
    data : FELLA_DATA
        Loaded database.
    approx : str
        Approximation: "simulation", "normality", "gamma", "t".
    t_df : float
        Degrees of freedom for t approximation.
    niter : int
        Number of permutations for simulation.

    Returns
    -------
    FELLA_USER
        Updated object.
    """
    if not is_FELLA_USER(object):
        print("'object' is not a FELLA_USER object. Returning original...")
        return object
    if not is_FELLA_DATA(data):
        print("'data' is not a FELLA_DATA object. Returning original...")
        return object
    if getStatus(data) != "loaded":
        print("'data' points to an empty FELLA_DATA object! Returning original...")
        return object

    print("Running diffusion...")
    comp_input = getInput(object)
    n_input = len(comp_input)
    if n_input == 0:
        print("Diffusion failed because there are no compounds in the input.")
        return object

    if approx == "simulation":
        print("Estimating p-values by simulation.")
        if len(getBackground(object)) == 0:
            comp_background = getCom(data, "compound")
        else:
            comp_background = getBackground(object)

        diffusion_matrix = getMatrix(data, "diffusion")
        if diffusion_matrix is None or diffusion_matrix.size <= 1:
            print("Diffusion matrix not loaded. Simulations will be slower...")
            graph = getGraph(data).as_undirected()
            n_nodes = graph.vcount()
            minusKI = np.array(graph.laplacian(normalized=False))
            id_pathway = list(data.keggdata.id.get("pathway", {}).values())
            diag = np.diag(minusKI).copy()
            diag[id_pathway] += 1.0
            np.fill_diagonal(minusKI, diag)

            generation = np.zeros(n_nodes)
            names = graph.vs["name"]
            idx_input = [names.index(c) for c in comp_input if c in names]
            generation[idx_input] = 1.0
            current_temp = la.solve(minusKI, generation)
            generation[idx_input] = 0.0

            null_temp = np.zeros((n_nodes, niter))
            bg_idx = [names.index(c) for c in comp_background if c in names]
            for dummy in range(niter):
                if (dummy + 1) % max(1, round(0.1 * niter)) == 0:
                    print(f"{round((dummy + 1) * 100 / niter)}%")
                sel = np.random.choice(bg_idx, size=n_input, replace=False)
                generation[sel] = 1.0
                null_temp[:, dummy] = la.solve(minusKI, generation)
                generation[sel] = 0.0

            pscores = np.zeros(n_nodes)
            for row in range(n_nodes):
                ecdf = stats.ecdf(null_temp[row, :])
                pscores[row] = (
                    (1 - ecdf.cdf.evaluate(current_temp[row])) * n_nodes + 1
                ) / (n_nodes + 1)
            pscores = dict(zip(names, pscores))
        else:
            names = getGraph(data).vs["name"]
            col_idx = [names.index(c) for c in comp_input if c in names]
            current_temp = diffusion_matrix[:, col_idx].sum(axis=1)
            bg_idx = [names.index(c) for c in comp_background if c in names]
            n_nodes = len(names)
            null_temp = np.zeros((n_nodes, niter))
            for dummy in range(niter):
                if (dummy + 1) % max(1, round(0.1 * niter)) == 0:
                    print(f"{round((dummy + 1) * 100 / niter)}%")
                sel = np.random.choice(bg_idx, size=n_input, replace=False)
                null_temp[:, dummy] = diffusion_matrix[:, sel].sum(axis=1)
            pscores = np.zeros(n_nodes)
            for row in range(n_nodes):
                ecdf = stats.ecdf(null_temp[row, :])
                pscores[row] = (
                    (1 - ecdf.cdf.evaluate(current_temp[row])) * n_nodes + 1
                ) / (n_nodes + 1)
            pscores = dict(zip(names, pscores))

    elif approx in ("normality", "gamma", "t"):
        print("Computing p-scores through the specified distribution.")
        if len(getBackground(object)) > 0:
            diffusion_matrix = getMatrix(data, "diffusion")
            if diffusion_matrix is None or diffusion_matrix.size <= 1:
                print(
                    "Diffusion matrix not loaded. Normality is not available "
                    "for custom background without it."
                )
                return object
            bg_names = getBackground(object)
            names = getGraph(data).vs["name"]
            col_idx = [names.index(c) for c in bg_names if c in names]
            background_matrix = diffusion_matrix[:, col_idx]
            RowSums = background_matrix.sum(axis=1)
            squaredRowSums = np.square(background_matrix).sum(axis=1)
            n_comp = background_matrix.shape[1]
        else:
            n_comp = len(getCom(data, "compound"))
            RowSums = getSums(data, "diffusion", squared=False)
            if RowSums is None or len(RowSums) == 0:
                print("RowSums not available. The normal approximation cannot be done.")
                return object
            squaredRowSums = getSums(data, "diffusion", squared=True)
            if squaredRowSums is None or len(squaredRowSums) == 0:
                print("squaredRowSums not available. The normal approximation cannot be done.")
                return object
            # Ensure they are numpy arrays
            if hasattr(RowSums, "values"):
                RowSums = RowSums.values
            if hasattr(squaredRowSums, "values"):
                squaredRowSums = squaredRowSums.values

        graph = getGraph(data).as_undirected()
        minusKI = np.array(graph.laplacian(normalized=False))
        id_pathway = list(data.keggdata.id.get("pathway", {}).values())
        diag = np.diag(minusKI).copy()
        diag[id_pathway] += 1.0
        np.fill_diagonal(minusKI, diag)

        names = graph.vs["name"]
        generation = np.zeros(len(names))
        idx_input = [names.index(c) for c in comp_input if c in names]
        generation[idx_input] = 1.0
        current_temp = la.solve(minusKI, generation)

        temp_means = RowSums * n_input / n_comp
        temp_vars = (
            n_input * (n_comp - n_input)
            / (n_comp * (n_comp - 1))
            * (squaredRowSums - (RowSums ** 2) / n_comp)
        )

        if approx == "normality":
            pscores = stats.norm.sf(current_temp, loc=temp_means, scale=np.sqrt(temp_vars))
        elif approx == "gamma":
            shape = temp_means ** 2 / temp_vars
            scale = temp_vars / temp_means
            pscores = stats.gamma.sf(current_temp, a=shape, scale=scale)
        elif approx == "t":
            pscores = stats.t.sf(
                (current_temp - temp_means) / np.sqrt(temp_vars),
                df=t_df,
            )
        pscores = dict(zip(names, pscores))
    else:
        raise ValueError(f"Unknown approx: {approx}")

    object.diffusion.pscores = pscores
    object.diffusion.approx = approx
    object.diffusion.niter = niter
    object.diffusion.valid = True
    print("Done.")
    return object
