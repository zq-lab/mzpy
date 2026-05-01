"""
Python port of FELLA::generateResultsTable
"""

from typing import Optional

import pandas as pd

from .AllClasses import FELLA_USER, FELLA_DATA
from .is_ import is_FELLA_DATA
from .get_ import getStatus, getPscores, getValid, getGraph, getInput
from .list_ import listCategories
from .AllArguments import checkArguments


def generateResultsTable(
    method: str = "diffusion",
    threshold: float = 0.05,
    plimit: int = 15,
    nlimit: int = 250,
    LabelLengthAtPlot: int = 45,
    capPscores: float = 1e-6,
    object: Optional[FELLA_USER] = None,
    data: Optional[FELLA_DATA] = None,
    **kwargs,
):
    """
    Generate a results table from an enrichment analysis.

    Returns
    -------
    pd.DataFrame or None
    """
    if not is_FELLA_DATA(data):
        raise ValueError("'data' is not a FELLA_DATA object")
    if getStatus(data) != "loaded":
        raise ValueError("'data' points to an empty FELLA_DATA object")

    check = checkArguments(
        method=method, threshold=threshold, plimit=plimit,
        nlimit=nlimit, LabelLengthAtPlot=LabelLengthAtPlot,
        object=object, data=data,
    )
    if not check["valid"]:
        raise ValueError("Bad argument when calling function 'generateResultsTable'.")

    valid = getValid(object, method)
    if valid is None or not valid:
        print(f"Method {method} has not been executed yet. Returning NULL...")
        return None

    if method == "hypergeom":
        print("Writing hypergeom results...")
        pscores = getPscores(object, "hypergeom")
        if not isinstance(pscores, pd.Series):
            pscores = pd.Series(pscores)
        pscores = pscores.sort_values()
        if pscores.iloc[0] >= threshold:
            print("No pathway is below the p-value threshold.")
            return None
        pscores = pscores[pscores < threshold].head(plimit)
        paths = pscores.index.tolist()
        pathnames = [str(data.keggdata.id2name.get(p, [""])[0]) for p in paths]
        pathhits = object.hypergeom.pathhits[paths]
        pathbackground = object.hypergeom.pathbackground[paths]
        out = pd.DataFrame({
            "KEGG id": paths,
            "KEGG name": pathnames,
            "CompoundHits": pathhits.values if hasattr(pathhits, "values") else list(pathhits),
            "CompoundsInPathway": pathbackground.values if hasattr(pathbackground, "values") else list(pathbackground),
            "p.value": pscores.values,
        })
        print("Done.")
        return out
    else:
        print(f"Writing {method} results...")
        pscores = getPscores(object, method)
        if not isinstance(pscores, pd.Series):
            pscores = pd.Series(pscores)
        pscores = pscores.sort_values()
        pscores[pscores < capPscores] = capPscores
        if pscores.iloc[0] >= threshold:
            print("No node is below the p-score threshold.")
            return None
        pscores = pscores[pscores < threshold].head(nlimit)
        pscores = pscores.sort_index()
        nodeIds = pscores.index.tolist()
        nodeNames = []
        for nid in nodeIds:
            ans = data.keggdata.id2name.get(nid, [""])
            if not ans:
                ans = ""
            else:
                ans = ans[0]
            if len(ans) > LabelLengthAtPlot:
                ans = ans[:LabelLengthAtPlot] + "..."
            nodeNames.append(ans)
        g = getGraph(data)
        nodeCom = [g.vs.find(name=nid)["com"] for nid in nodeIds]
        nodeTypes = [listCategories()[c - 1] for c in nodeCom]
        out = pd.DataFrame({
            "KEGG id": nodeIds,
            "Entry type": nodeTypes,
            "KEGG name": nodeNames,
            "p.score": pscores.values,
        })
        out = out.sort_values(by="Entry type")
        out.reset_index(drop=True, inplace=True)
        print("Done.")
        return out
