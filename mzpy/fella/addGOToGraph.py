"""
Python port of FELLA::addGOToGraph

Uses the BioMart REST API (via bioservices or raw requests) to map
entrez genes to GO terms, and goatools (when available) to compute
semantic similarities.  If the GO OBO file cannot be fetched, only
the GO labels are appended.
"""

import os
from typing import Optional, List, Dict

import igraph as ig


# ---------------------------------------------------------------------------
# BioMart helper
# ---------------------------------------------------------------------------

def _query_biomart(dataset: str, attributes: List[str], filters: Dict[str, str]) -> List[str]:
    """Query BioMart via raw POST (bioservices wrapper is too variable)."""
    try:
        from bioservices import BioMart
        bm = BioMart()
    except Exception:
        bm = None

    xml = (
        '<?xml version="1.0" encoding="UTF-8"?>\n'
        '<!DOCTYPE Query>\n'
        '<Query virtualSchemaName="default" formatter="TSV" header="0" '
        'uniqueRows="0" count="" datasetConfigVersion="0.6">\n'
        f'  <Dataset name="{dataset}" interface="default">\n'
    )
    for fname, fvalue in filters.items():
        # BioMart expects comma-separated values in the value attribute
        val = fvalue if isinstance(fvalue, str) else ",".join(str(x) for x in fvalue)
        xml += f'    <Filter name="{fname}" value="{val}"/>\n'
    for attr in attributes:
        xml += f'    <Attribute name="{attr}"/>\n'
    xml += (
        '  </Dataset>\n'
        '</Query>\n'
    )

    if bm is not None:
        try:
            res = bm.query(xml)
            return [line for line in res.strip().splitlines() if line]
        except Exception:
            pass

    # fallback to requests
    import requests
    url = "http://www.ensembl.org/biomart/martservice"
    resp = requests.post(url, data={"query": xml}, timeout=60)
    resp.raise_for_status()
    return [line for line in resp.text.strip().splitlines() if line]


# ---------------------------------------------------------------------------
# GO semantic similarity helper (goatools)
# ---------------------------------------------------------------------------

def _try_load_godag(obo_path: str = "go-basic.obo"):
    """Try to load a local GODag; download if missing and network allows."""
    try:
        from goatools.obo_parser import GODag
    except Exception:
        return None
    if not os.path.exists(obo_path):
        try:
            from goatools.base import download_go_basic_obo
            download_go_basic_obo(obo_path)
        except Exception as exc:
            print(
                f"Could not download {obo_path} ({exc}). "
                "GO semantic similarity will be skipped."
            )
            return None
    try:
        return GODag(obo_path)
    except Exception:
        return None


def _max_resnik_sim(query_go: str, go_list: List[str], godag, termcounts) -> Optional[float]:
    """Return the maximum Resnik similarity between query_go and any term in go_list."""
    from goatools.semantic import resnik_sim
    max_sim = 0.0
    for g in go_list:
        try:
            sim = resnik_sim(query_go, g, godag, termcounts)
            if sim is not None and sim > max_sim:
                max_sim = sim
        except Exception:
            pass
    return max_sim if max_sim > 0 else None


# ---------------------------------------------------------------------------
# Main function
# ---------------------------------------------------------------------------

def addGOToGraph(
    graph: Optional[ig.Graph] = None,
    GOterm: Optional[str] = None,
    godata_options: Optional[dict] = None,
    mart_options: Optional[dict] = None,
):
    """
    Add GO labels (and optionally semantic similarity) to a graph.

    Parameters
    ----------
    graph : igraph.Graph
        Input graph (must contain vertex attribute 'entrez').
    GOterm : str or None
        Query GO term, e.g. "GO:0005739".  If None, only GO labels
        are appended without similarity values.
    godata_options : dict or None
        Ignored in the Python port (kept for API compatibility).
    mart_options : dict or None
        BioMart dataset specification.  Expected keys:
        ``dataset`` (default "hsapiens_gene_ensembl").

    Returns
    -------
    igraph.Graph
        The graph with GO attributes added.
    """
    if graph is None or not isinstance(graph, ig.Graph):
        print("'graph' is not an igraph object. Leaving it as it was...")
        return graph
    if graph.vcount() == 0:
        return graph

    # Default mart options
    if mart_options is None:
        mart_options = {"dataset": "hsapiens_gene_ensembl"}
    dataset = mart_options.get("dataset", "hsapiens_gene_ensembl")

    # Collect all entrez ids from the graph
    entrez_all = set()
    for v in graph.vs:
        genes = v.attributes().get("entrez", [])
        if genes:
            entrez_all.update(str(g) for g in genes)
    entrez_all = sorted(entrez_all)
    if not entrez_all:
        print("No entrez gene ids found in the graph. GO labels cannot be added.")
        return graph

    # Query BioMart
    print("Querying BioMart for GO annotations...")
    try:
        lines = _query_biomart(
            dataset=dataset,
            attributes=["entrezgene_id", "go_id"],
            filters={"entrezgene_id": ",".join(entrez_all)},
        )
    except Exception as exc:
        print(f"BioMart query failed ({exc}). GO labels will not be added.")
        return graph

    # Parse TSV
    entrez2go: Dict[str, List[str]] = {}
    for line in lines:
        parts = line.split("\t")
        if len(parts) < 2:
            continue
        eid, go_id = parts[0], parts[1].strip()
        if go_id.startswith("GO:"):
            entrez2go.setdefault(eid, []).append(go_id)

    # Attach GO lists to vertices
    for v in graph.vs:
        genes = v.attributes().get("entrez", [])
        gos = set()
        for g in genes:
            gos.update(entrez2go.get(str(g), []))
        v["GO"] = sorted(gos)

    if GOterm is None:
        print(
            "Null GOterm provided to addGOToGraph. "
            "Only the GO labels will be added. "
            "To include similarity values as well, please specify a GOterm."
        )
        return graph

    # Compute semantic similarities
    print("Computing GO semantic similarities...")
    godag = _try_load_godag()
    if godag is None:
        print(
            "GO OBO file not available. "
            "Similarity values will not be computed."
        )
        for v in graph.vs:
            v["GO.simil"] = None
        return graph

    # Build TermCounts from the graph's GO annotations (approximation)
    all_go_lists = [v["GO"] for v in graph.vs if v["GO"]]
    if not all_go_lists:
        print("No GO annotations found. Similarity values will not be computed.")
        for v in graph.vs:
            v["GO.simil"] = None
        return graph

    try:
        from goatools.semantic import TermCounts
        assocs = {i: g for i, g in enumerate(all_go_lists)}
        termcounts = TermCounts(godag, assocs)
    except Exception as exc:
        print(f"Could not build TermCounts ({exc}). Similarity values skipped.")
        for v in graph.vs:
            v["GO.simil"] = None
        return graph

    for v in graph.vs:
        gos = v["GO"]
        if gos:
            sim = _max_resnik_sim(GOterm, gos, godag, termcounts)
            v["GO.simil"] = sim
        else:
            v["GO.simil"] = None

    print("Done.")
    return graph
