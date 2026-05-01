"""
Python port of FELLA::get* extractor functions.
"""

from typing import Any, List, Optional

from .AllClasses import FELLA_DATA, FELLA_USER


def getBackground(object: FELLA_USER) -> List[str]:
    """Return the compounds defined as background."""
    return object.userinput.metabolitesbackground


def getCom(data: FELLA_DATA, level, format: str = "name"):
    """
    Return all nodes from a level/community of the KEGG graph.
    format="name" returns KEGG IDs; format="id" returns vertex indices.
    """
    # level can be int 1..5 or str like "pathway"
    if isinstance(level, int):
        cats = ["pathway", "module", "enzyme", "reaction", "compound"]
        level = cats[level - 1]
    if format == "name":
        return list(data.keggdata.id.get(level, {}).keys())
    if format == "id":
        return data.keggdata.id.get(level, {})
    raise ValueError("'format' must be 'name' or 'id'.")


def getExcluded(object: FELLA_USER) -> List[str]:
    """Return the compounds in the input that were not mapped."""
    return object.userinput.excluded


def getGraph(data: FELLA_DATA):
    """Return the KEGG graph (igraph object)."""
    return data.keggdata.graph


def getInfo(data: FELLA_DATA) -> str:
    """Return the KEGG version info used to build the database."""
    g = getGraph(data)
    return g.get("comment", "") if g else ""


def getInput(object: FELLA_USER) -> List[str]:
    """Return the metabolites specified by the user."""
    return object.userinput.metabolites


def getMatrix(data: FELLA_DATA, method: str):
    """Return the matrix for the desired methodology."""
    if method == "hypergeom":
        return data.hypergeom.matrix
    if method == "diffusion":
        return data.diffusion.matrix
    if method == "pagerank":
        return data.pagerank.matrix
    raise ValueError("method must be one of 'hypergeom', 'diffusion', 'pagerank'.")


def getName(data: FELLA_DATA, id: List[str]):
    """Map KEGG identifiers to KEGG names."""
    out = {}
    for i in id:
        out[i] = data.keggdata.id2name.get(i, [""])
    return out


def getPscores(object: FELLA_USER, method: str):
    """Return p-scores / p-values from the desired methodology."""
    if method == "hypergeom":
        return object.hypergeom.pvalues
    if method == "diffusion":
        return object.diffusion.pscores
    if method == "pagerank":
        return object.pagerank.pscores
    raise ValueError("method must be one of 'hypergeom', 'diffusion', 'pagerank'.")


def getPvaluesSize(data: FELLA_DATA):
    """Return the matrix containing p-values by CC size."""
    return data.keggdata.pvalues_size


def getSums(data: FELLA_DATA, method: str, squared: bool = False):
    """Return rowSums or squaredRowSums for a method."""
    if method == "diffusion":
        if squared:
            return data.diffusion.squaredRowSums
        return data.diffusion.rowSums
    if method == "pagerank":
        if squared:
            return data.pagerank.squaredRowSums
        return data.pagerank.rowSums
    raise ValueError("method must be 'diffusion' or 'pagerank'.")


def getValid(object: FELLA_USER, method: str) -> Optional[bool]:
    """Return the 'valid' slot for a method."""
    if method == "hypergeom":
        return object.hypergeom.valid
    if method == "diffusion":
        return object.diffusion.valid
    if method == "pagerank":
        return object.pagerank.valid
    raise ValueError("method must be one of 'hypergeom', 'diffusion', 'pagerank'.")


def getStatus(data: FELLA_DATA) -> str:
    """Return the 'status' slot of the KEGG data."""
    return data.keggdata.status
