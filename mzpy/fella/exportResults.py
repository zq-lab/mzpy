"""
Python port of FELLA::exportResults
"""

from typing import Optional

from .AllClasses import FELLA_USER, FELLA_DATA
from .is_ import is_FELLA_DATA
from .get_ import getStatus
from .generateResultsTable import generateResultsTable
from .generateEnzymesTable import generateEnzymesTable
from .generateResultsGraph import generateResultsGraph


def exportResults(
    format: str = "csv",
    file: str = "myOutput",
    method: str = "diffusion",
    object: Optional[FELLA_USER] = None,
    data: Optional[FELLA_DATA] = None,
    **kwargs,
):
    """
    Export enrichment results to a file.

    Parameters
    ----------
    format : str
        One of "csv", "enzyme", or any igraph-supported graph format.
    file : str
        Output file path.
    method : str
        Enrichment method.
    object : FELLA_USER
        Results object.
    data : FELLA_DATA
        Database.
    """
    if not is_FELLA_DATA(data):
        raise ValueError("'data' is not a FELLA_DATA object")
    if getStatus(data) != "loaded":
        raise ValueError("'data' points to an empty FELLA_DATA object")

    if format == "csv":
        print("Exporting to a csv file...")
        df = generateResultsTable(
            method=method, object=object, data=data, **kwargs
        )
        if df is not None:
            df.to_csv(file, index=False)
    elif format == "enzyme":
        if method in ("diffusion", "pagerank"):
            df = generateEnzymesTable(
                method=method, object=object, data=data, **kwargs
            )
            if df is not None:
                df.to_csv(file, index=False)
        else:
            raise ValueError(
                f"Enzymes are only reported in diffusion and pagerank, "
                f"but not in method {method}"
            )
    else:
        graph = generateResultsGraph(
            method=method, object=object, data=data, **kwargs
        )
        if graph is None:
            print("Graph is empty, nothing to export.")
            return
        if format == "igraph":
            import pickle
            with open(file, "wb") as f:
                pickle.dump({"graph": graph}, f)
            print("Exported to pickle file using igraph object...")
        else:
            graph.write(file, format=format)
            print(f"Exporting to the format {format} using igraph...")

    print("Done")
