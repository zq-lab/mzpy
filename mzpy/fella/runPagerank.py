"""
Python port of FELLA::runPagerank
"""

from typing import Optional

import numpy as np
from scipy import stats

from .AllClasses import FELLA_USER, FELLA_DATA
from .is_ import is_FELLA_USER, is_FELLA_DATA
from .get_ import (
    getStatus, getInput, getBackground, getMatrix, getCom,
    getGraph, getSums,
)


def runPagerank(
    object: Optional[FELLA_USER] = None,
    data: Optional[FELLA_DATA] = None,
    approx: str = "normality",
    dampingFactor: float = 0.85,
    t_df: float = 10.0,
    niter: int = 1000,
) -> FELLA_USER:
    """
    Perform PageRank-based enrichment.

    Parameters
    ----------
    object : FELLA_USER
        User object with mapped metabolites.
    data : FELLA_DATA
        Loaded database.
    approx : str
        Approximation: "simulation", "normality", "gamma", "t".
    dampingFactor : float
        Damping factor for PageRank.
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
        print("'data' points to an empty FELLA.DATA object! Returning original...")
        return object

    print("Running PageRank...")
    d = dampingFactor
    comp_input = getInput(object)
    n_input = len(comp_input)
    if n_input == 0:
        print("PageRank failed because there are no compounds in the input.")
        return object

    if approx == "simulation":
        print("Estimating p-values by simulation.")
        if len(getBackground(object)) == 0:
            comp_background = getCom(data, "compound")
        else:
            comp_background = getBackground(object)

        pagerank_matrix = getMatrix(data, "pagerank")
        if pagerank_matrix is None or pagerank_matrix.size <= 1:
            print("PageRank matrix not loaded. Simulations may be a bit slower...")
            graph = getGraph(data)
            n_nodes = graph.vcount()
            names = graph.vs["name"]
            prior = np.zeros(n_nodes)
            idx_input = [names.index(c) for c in comp_input if c in names]
            prior[idx_input] = 1.0
            # igraph page_rank personalised uses restart vector
            current_score = np.array(
                graph.personalized_pagerank(damping=d, reset_vertices=prior)
            )
            prior[idx_input] = 0.0

            bg_idx = [names.index(c) for c in comp_background if c in names]
            null_score = np.zeros((n_nodes, niter))
            for dummy in range(niter):
                if (dummy + 1) % max(1, round(0.1 * niter)) == 0:
                    print(f"{round((dummy + 1) * 100 / niter)}%")
                sel = np.random.choice(bg_idx, size=n_input, replace=False)
                prior[sel] = 1.0
                null_score[:, dummy] = np.array(
                    graph.personalized_pagerank(damping=d, reset_vertices=prior)
                )
                prior[sel] = 0.0

            pscores = np.zeros(n_nodes)
            for row in range(n_nodes):
                ecdf = stats.ecdf(null_score[row, :])
                pscores[row] = (
                    (1 - ecdf.cdf.evaluate(current_score[row])) * n_nodes + 1
                ) / (n_nodes + 1)
            pscores = dict(zip(names, pscores))
        else:
            names = getGraph(data).vs["name"]
            col_idx = [names.index(c) for c in comp_input if c in names]
            current_score = pagerank_matrix[:, col_idx].sum(axis=1)
            bg_idx = [names.index(c) for c in comp_background if c in names]
            n_nodes = len(names)
            null_score = np.zeros((n_nodes, niter))
            for dummy in range(niter):
                if (dummy + 1) % max(1, round(0.1 * niter)) == 0:
                    print(f"{round((dummy + 1) * 100 / niter)}%")
                sel = np.random.choice(bg_idx, size=n_input, replace=False)
                null_score[:, dummy] = pagerank_matrix[:, sel].sum(axis=1)
            pscores = np.zeros(n_nodes)
            for row in range(n_nodes):
                ecdf = stats.ecdf(null_score[row, :])
                pscores[row] = (
                    (1 - ecdf.cdf.evaluate(current_score[row])) * n_nodes + 1
                ) / (n_nodes + 1)
            pscores = dict(zip(names, pscores))

    elif approx in ("normality", "gamma", "t"):
        print("Computing p-scores through the specified distribution.")
        if len(getBackground(object)) > 0:
            pagerank_matrix = getMatrix(data, "pagerank")
            if pagerank_matrix is None or pagerank_matrix.size <= 1:
                print(
                    "PageRank matrix not loaded. Normality is not available "
                    "for custom background without it."
                )
                return object
            bg_names = getBackground(object)
            names = getGraph(data).vs["name"]
            col_idx = [names.index(c) for c in bg_names if c in names]
            background_matrix = pagerank_matrix[:, col_idx]
            RowSums = background_matrix.sum(axis=1)
            squaredRowSums = np.square(background_matrix).sum(axis=1)
            n_comp = background_matrix.shape[1]
        else:
            n_comp = len(getCom(data, "compound"))
            RowSums = getSums(data, "pagerank", squared=False)
            if RowSums is None or len(RowSums) == 0:
                print("RowSums not available. The normal approximation cannot be done.")
                return object
            squaredRowSums = getSums(data, "pagerank", squared=True)
            if squaredRowSums is None or len(squaredRowSums) == 0:
                print("squaredRowSums not available. The normal approximation cannot be done.")
                return object
            if hasattr(RowSums, "values"):
                RowSums = RowSums.values
            if hasattr(squaredRowSums, "values"):
                squaredRowSums = squaredRowSums.values

        graph = getGraph(data)
        names = graph.vs["name"]
        prior = np.zeros(len(names))
        idx_input = [names.index(c) for c in comp_input if c in names]
        prior[idx_input] = 1.0
        current_score = np.array(
            graph.personalized_pagerank(damping=d, reset_vertices=prior)
        )

        score_means = RowSums / n_comp
        score_vars = (
            (n_comp - n_input)
            / (n_input * n_comp * (n_comp - 1))
            * (squaredRowSums - (RowSums ** 2) / n_comp)
        )

        if approx == "normality":
            pscores = stats.norm.sf(current_score, loc=score_means, scale=np.sqrt(score_vars))
        elif approx == "gamma":
            shape = score_means ** 2 / score_vars
            scale = score_vars / score_means
            pscores = stats.gamma.sf(current_score, a=shape, scale=scale)
        elif approx == "t":
            pscores = stats.t.sf(
                (current_score - score_means) / np.sqrt(score_vars),
                df=t_df,
            )
        pscores = dict(zip(names, pscores))
    else:
        raise ValueError(f"Unknown approx: {approx}")

    object.pagerank.pscores = pscores
    object.pagerank.approx = approx
    object.pagerank.niter = niter
    object.pagerank.valid = True
    print("Done.")
    return object
