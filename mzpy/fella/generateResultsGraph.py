"""
Python port of FELLA::generateResultsGraph
"""

from typing import Optional

import numpy as np
import igraph as ig

from .AllClasses import FELLA_USER, FELLA_DATA
from .is_ import is_FELLA_DATA
from .get_ import getStatus, getPscores, getValid, getGraph, getInput, getMatrix, getPvaluesSize
from .AllArguments import checkArguments


def generateResultsGraph(
    method: str = "diffusion",
    threshold: float = 0.05,
    plimit: int = 15,
    nlimit: int = 250,
    thresholdConnectedComponent: float = 0.05,
    LabelLengthAtPlot: int = 22,
    minDegree: int = 0,
    object: Optional[FELLA_USER] = None,
    data: Optional[FELLA_DATA] = None,
    **kwargs,
):
    """
    Generate a sub-network with the nodes with the lowest p-score.

    Parameters
    ----------
    method : str
        One of "hypergeom", "diffusion", "pagerank".
    threshold : float
        p-score threshold.
    plimit : int
        Max pathways for hypergeom.
    nlimit : int
        Max nodes for diffusion/pagerank.
    thresholdConnectedComponent : float
        Filter tiny connected components.
    LabelLengthAtPlot : int
        Max label length.
    object : FELLA_USER
        Enrichment results.
    data : FELLA_DATA
        Database.

    Returns
    -------
    igraph.Graph or None
        Sub-network graph.
    """
    if not is_FELLA_DATA(data):
        raise ValueError("'data' is not a FELLA_DATA object")
    if getStatus(data) != "loaded":
        raise ValueError("'data' points to an empty FELLA_DATA object")

    check = checkArguments(
        method=method,
        threshold=threshold,
        plimit=plimit,
        nlimit=nlimit,
        thresholdConnectedComponent=thresholdConnectedComponent,
        LabelLengthAtPlot=LabelLengthAtPlot,
        object=object,
        data=data,
    )
    if not check["valid"]:
        raise ValueError("Bad argument when calling function 'generateResultsGraph'.")

    valid = getValid(object, method)
    if valid is None or not valid:
        print(f"Method {method} has not been executed yet. Returning NULL...")
        return None

    if method == "hypergeom":
        pscores = getPscores(object, "hypergeom")
        g_data = getGraph(data)
        n_paths = sum(1 for p in pscores.values() if p < threshold)
        if n_paths < 1:
            print("Graph is empty. None of the pathways is significant.")
            return None
        else:
            sorted_pscores = sorted(pscores.items(), key=lambda x: x[1])
            path_hypergeom = [k for k, v in sorted_pscores if v < threshold][:plimit]

        comp_hypergeom = [c for c in getInput(object) if c in data.hypergeom.row_names]

        incidence = getMatrix(data, "hypergeom")[
            [data.hypergeom.row_names.index(c) for c in comp_hypergeom if c in data.hypergeom.row_names],
            :,
        ][:, [data.hypergeom.col_names.index(p) for p in path_hypergeom if p in data.hypergeom.col_names]]

        # Build bipartite graph
        # incidence is compound x pathway
        # igraph.Graph.Incidence expects incidence matrix (row=type1, col=type2)
        graph_bipartite = ig.Graph.Incidence(incidence.toarray().tolist(), directed=False)
        # Add names
        names = comp_hypergeom + path_hypergeom
        for i, v in enumerate(graph_bipartite.vs):
            v["name"] = names[i]
        # delete isolated
        graph_bipartite = graph_bipartite.subgraph([v.index for v in graph_bipartite.vs if graph_bipartite.degree(v.index) > 0])

        # filter by minimum degree
        if minDegree > 0:
            graph_bipartite = graph_bipartite.subgraph(
                [v.index for v in graph_bipartite.vs if graph_bipartite.degree(v.index) >= minDegree]
            )

        # com attribute
        name_to_com = {v["name"]: v["com"] for v in g_data.vs}
        for v in graph_bipartite.vs:
            v["com"] = name_to_com.get(v["name"], None)
        # label
        name_to_name = {v["name"]: v["NAME"] for v in g_data.vs}
        for v in graph_bipartite.vs:
            raw = name_to_name.get(v["name"], [""])
            if not raw:
                raw = [""]
            first = raw[0]
            v["label"] = first[:LabelLengthAtPlot] + ("..." if len(first) > LabelLengthAtPlot else "")
        return graph_bipartite

    else:
        pscores = getPscores(object, method)
        inp = getInput(object)
        nodes = [k for k, v in sorted(pscores.items(), key=lambda x: x[1]) if v < threshold]
        n_nodes = len(nodes)
        if n_nodes < 1:
            print("Graph is empty.")
            return None
        elif n_nodes > nlimit:
            print(f"{n_nodes} nodes below the threshold have been limited to {nlimit} nodes.")
            nodes = nodes[:nlimit]

        graph = getGraph(data).induced_subgraph(nodes)
        if method == "diffusion":
            graph = graph.as_undirected()

        # labels
        for v in graph.vs:
            raw = v["NAME"] if "NAME" in v.attributes() else [""]
            if not raw:
                raw = [""]
            first = raw[0]
            v["label"] = first[:LabelLengthAtPlot] + ("..." if len(first) > LabelLengthAtPlot else "")
            v["input"] = v["name"] in inp

        # Filter small CCs
        n_nodes_graph = graph.vcount()
        graph_clust = graph.clusters(mode="weak")

        if n_nodes_graph > 250:
            print(
                f"The number of nodes of the whole solution ({n_nodes_graph}) exceeds 250. "
                f"Small connected components will be filtered using 250 nodes instead."
            )
            n_nodes_graph = 250

        pvalues_size = getPvaluesSize(data)
        if pvalues_size is not None and pvalues_size.shape[1] >= n_nodes_graph:
            tab_significant = pvalues_size[:, n_nodes_graph - 1]
            csize_significant = next(
                (i for i, p in enumerate(tab_significant) if p <= thresholdConnectedComponent),
                None,
            )
            if csize_significant is None:
                print(
                    "None of the connected components are below the "
                    "'thresholdConnectedComponent'. Returning empty graph..."
                )
                return ig.Graph()
            csize_significant += 1  # 1-based
            clust_select = [i for i, sz in enumerate(graph_clust.sizes()) if sz >= csize_significant]
            if not clust_select:
                print(
                    "No connected component meets the size threshold "
                    f"(>={csize_significant}). Returning empty graph..."
                )
                return ig.Graph()
            graph = graph.induced_subgraph(
                [v for v in graph.vs if graph_clust.membership[v.index] in clust_select]
            )

        # filter by minimum degree
        if minDegree > 0:
            graph = graph.subgraph(
                [v.index for v in graph.vs if graph.degree(v.index) >= minDegree]
            )

        return graph
