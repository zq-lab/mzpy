"""
Python port of FELLA::generateEnzymesTable
"""

from typing import Optional

import pandas as pd

from .AllClasses import FELLA_USER, FELLA_DATA
from .is_ import is_FELLA_DATA
from .get_ import getStatus, getPscores, getValid, getGraph, getCom, getName
from .AllArguments import checkArguments


def generateEnzymesTable(
    method: str = "diffusion",
    threshold: float = 0.05,
    nlimit: int = 250,
    LabelLengthAtPlot: int = 45,
    capPscores: float = 1e-6,
    mart_options: Optional[dict] = None,
    object: Optional[FELLA_USER] = None,
    data: Optional[FELLA_DATA] = None,
    **kwargs,
):
    """
    Generate a table with enzyme families and their annotated genes.

    Returns
    -------
    pd.DataFrame or None
    """
    check = checkArguments(
        method=method, threshold=threshold, nlimit=nlimit,
        LabelLengthAtPlot=LabelLengthAtPlot,
        object=object, data=data,
    )
    if not check["valid"]:
        raise ValueError("Bad argument when calling function 'generateEnzymesTable'.")

    if getStatus(data) != "loaded":
        raise ValueError("'data' points to an empty FELLA_DATA object")

    if method not in ("diffusion", "pagerank"):
        print(f"Method should be 'diffusion' or 'pagerank', but it is {method}. Returning NULL...")
        return None

    valid = getValid(object, method)
    if valid is None or not valid:
        print(f"Method {method} has not been executed yet. Returning NULL...")
        return None

    print(f"Writing {method} enzymes...")
    id_ec = getCom(data, level=3, format="id")
    pscores_ec = {k: v for k, v in getPscores(object, method).items() if k in id_ec}
    if not pscores_ec:
        print("No enzyme is below the p-score threshold.")
        return None
    pscores_ec = pd.Series(pscores_ec).sort_values()
    pscores_ec[pscores_ec < capPscores] = capPscores
    if pscores_ec.iloc[0] >= threshold:
        print("No enzyme is below the p-score threshold.")
        return None
    nodePscores = pscores_ec[pscores_ec < threshold].head(nlimit)
    nodeIds = nodePscores.index.tolist()

    nodeNames = []
    for nid in nodeIds:
        ans = getName(data, [nid]).get(nid, [""])
        if not ans:
            ans = ""
        else:
            ans = ans[0]
        if len(ans) > LabelLengthAtPlot:
            ans = ans[:LabelLengthAtPlot] + "..."
        nodeNames.append(ans)

    g = getGraph(data)
    nodeGenes = []
    for nid in nodeIds:
        v = g.vs.find(name=nid)
        genes = v.attributes().get("entrez", [])
        nodeGenes.append(";".join(str(g) for g in genes) if genes else "")

    out = pd.DataFrame({
        "EC_number": nodeIds,
        "p.score": nodePscores.values,
        "EC_name": nodeNames,
        "Genes": nodeGenes,
    })
    print("Done.")
    return out
