"""
Python port of FELLA::AllMethods (show and plot S4 methods).

Python does not have S4, so we provide free functions ``show`` and
``plot`` that mimic the original behaviour.
"""

from typing import Any, Optional

import pandas as pd

from .AllClasses import FELLA_DATA, FELLA_USER
from .get_ import getGraph, getMatrix, getSums, getValid, getPscores, getInput, getBackground, getStatus
from .is_ import is_FELLA_DATA, is_FELLA_USER
from .list_ import listCategories, listMethods
from .generateResultsGraph import generateResultsGraph
from .plotBipartite import plotBipartite
from .plotGraph import plotGraph
from .AllArguments import checkArguments


def show(obj: Any) -> Optional[pd.DataFrame]:
    """
    Print a summary of a FELLA_DATA or FELLA_USER object.

    Returns a small DataFrame with top hits if obj is FELLA_USER.
    """
    if is_FELLA_DATA(obj):
        _show_data(obj)
        return None
    elif is_FELLA_USER(obj):
        return _show_user(obj)
    else:
        print("Object is neither FELLA_DATA nor FELLA_USER.")
        return None


def _show_data(obj: FELLA_DATA):
    print("---------------------")
    print("General data:")
    g = getGraph(obj)
    if g is None or g.vcount() == 0:
        print("- KEGG graph not loaded.")
    else:
        print("- KEGG graph:")
        print(f"  * Nodes: {g.vcount()}")
        print(f"  * Edges: {g.ecount()}")
        # density
        if g.vcount() > 1:
            dens = g.ecount() / (g.vcount() * (g.vcount() - 1))
        else:
            dens = 0.0
        print(f"  * Density: {dens:.4f}")
        cats = listCategories()
        coms = [v["com"] for v in g.vs]
        tab = {cat: sum(1 for c in coms if c == i + 1) for i, cat in enumerate(cats)}
        for cat, cnt in tab.items():
            print(f"    + {cat} [{cnt}]")
    if len(obj.keggdata.id2name) == 0:
        print("- KEGG names not loaded.")
    else:
        print("- KEGG names are ready.")
    print("---------------------")

    print("Hypergeometric test:")
    mat = getMatrix(obj, "hypergeom")
    if mat is None or mat.size <= 1:
        print("- Matrix not loaded.")
    else:
        print(f"- Matrix is ready\n  * Dim: {mat.shape[0]}x{mat.shape[1]}")
    print("---------------------")

    print("Heat diffusion:")
    mat = getMatrix(obj, "diffusion")
    if mat is None or mat.size <= 1:
        print("- Matrix not loaded.")
    else:
        print(f"- Matrix is ready\n  * Dim: {mat.shape[0]}x{mat.shape[1]}")
    rs = getSums(obj, "diffusion", squared=False)
    if rs is None or len(rs) == 0:
        print("- RowSums not loaded.")
    else:
        print("- RowSums are ready.")
    print("---------------------")

    print("PageRank:")
    mat = getMatrix(obj, "pagerank")
    if mat is None or mat.size <= 1:
        print("- Matrix not loaded.")
    else:
        print(f"- Matrix is ready\n  * Dim: {mat.shape[0]}x{mat.shape[1]}")
    rs = getSums(obj, "pagerank", squared=False)
    if rs is None or len(rs) == 0:
        print("- RowSums not loaded.")
    else:
        print("- RowSums are ready.")
    print("---------------------")


def _show_user(obj: FELLA_USER) -> dict:
    print("---------------------")
    inp = getInput(obj)
    print(f"Compounds in the input: {len(inp)}")
    if len(inp) > 0:
        print(inp)
    bg = getBackground(obj)
    if len(bg) == 0:
        print("Background compounds: all available compounds (default)")
    else:
        print(f"Background compounds: {len(bg)}")
    print("---------------------")

    for method in listMethods():
        valid = getValid(obj, method)
        if valid is None:
            print(f"{method}: not performed")
        elif not valid:
            print(f"{method}: error during execution")
        else:
            print(f"{method}: ready.")
            if method != "hypergeom":
                pscores = getPscores(obj, method)
                if isinstance(pscores, dict):
                    below = sum(1 for v in pscores.values() if v < 0.05)
                else:
                    below = sum(1 for v in pscores if v < 0.05)
                print(f"  P-scores under 0.05: {below}")
        print("---------------------")

    # Return small data frames with top hits
    ans = {}
    for method in listMethods():
        valid = getValid(obj, method)
        if valid is None or not valid:
            ans[method] = "Not performed"
            continue
        if method == "hypergeom":
            pval = getPscores(obj, "hypergeom")
            if isinstance(pval, pd.Series):
                pval = pval.sort_values().head(30)
            else:
                pval = pd.Series(pval).sort_values().head(30)
            ans[method] = pd.DataFrame({
                "Description": pval.index,
                "p.value": pval.values,
            })
        else:
            pscores = getPscores(obj, method)
            if isinstance(pscores, dict):
                pscores = pd.Series(pscores)
            else:
                pscores = pd.Series(pscores)
            pscores = pscores.sort_values().head(50)
            pscores = pscores.sort_index()
            pscores_str = pscores.apply(lambda x: "<1e-6" if x < 1e-6 else f"{x:.4g}")
            ans[method] = pd.DataFrame({
                "KEGG id": pscores.index,
                "p.score": pscores_str.values,
            })
    return ans


def plot(
    x: FELLA_USER,
    method: str = "hypergeom",
    threshold: float = 0.05,
    plimit: int = 15,
    nlimit: int = 250,
    layout: bool = False,
    thresholdConnectedComponent: float = 0.05,
    LabelLengthAtPlot: int = 22,
    data: Optional[FELLA_DATA] = None,
    **kwargs,
):
    """
    Plot results from a FELLA_USER object.

    Parameters
    ----------
    x : FELLA_USER
        Enrichment results.
    method : str
        Method to plot.
    threshold : float
        p-score threshold.
    plimit : int
        Max pathways for hypergeom.
    nlimit : int
        Max nodes for diffusion/pagerank.
    layout : bool
        If True, return layout DataFrame.
    thresholdConnectedComponent : float
        CC size filter.
    LabelLengthAtPlot : int
        Max label length.
    data : FELLA_DATA
        Database.
    """
    check = checkArguments(
        method=method, threshold=threshold, plimit=plimit, nlimit=nlimit,
        layout=layout, thresholdConnectedComponent=thresholdConnectedComponent,
        LabelLengthAtPlot=LabelLengthAtPlot, object=x, data=data,
    )
    if not check["valid"]:
        raise ValueError("Bad argument when calling function 'plot' in FELLA.")

    if getStatus(data) != "loaded":
        raise ValueError("'data' points to an empty FELLA_DATA object")

    valid = getValid(x, method)
    if valid is None or not valid:
        raise ValueError(f"Results from {method} are not ready yet.")

    graph = generateResultsGraph(
        method=method, threshold=threshold, plimit=plimit, nlimit=nlimit,
        thresholdConnectedComponent=thresholdConnectedComponent,
        LabelLengthAtPlot=LabelLengthAtPlot, object=x, data=data, **kwargs
    )

    if graph is None:
        print("Graph is empty, nothing to plot.")
        return None

    if method == "hypergeom":
        return plotBipartite(graph=graph, layout=layout, **kwargs)
    else:
        return plotGraph(graph=graph, layout=layout, **kwargs)
