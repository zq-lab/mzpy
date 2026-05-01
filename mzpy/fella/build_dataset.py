"""
Build a FELLA dataset (KEGG graph + HDF5 database) for a given organism.

This module provides a single convenience function that wraps the two-step
process of fetching KEGG data and compiling the internal matrices.
"""

import os
import pickle
from typing import List, Optional

from .buildGraphFromKEGGREST import buildGraphFromKEGGREST
from .buildDataFromGraph import buildDataFromGraph


def build_dataset(
    organism: str,
    output_dir: str = ".",
    overwrite: bool = False,
    filter_path: Optional[List[str]] = None,
    matrices: Optional[List[str]] = None,
    normality: Optional[List[str]] = None,
    damping_factor: float = 0.85,
    niter: int = 100,
) -> dict:
    """
    Build and save a FELLA dataset for the specified organism.

    The function performs two main steps:
    1. Query the KEGG REST API to build a curated igraph object.
    2. Compute internal matrices (hypergeometric, diffusion, pagerank) and
       store everything in an HDF5 file.

    Before writing, existing files are checked.  If a file already exists and
    *overwrite* is ``False``, a ``FileExistsError`` is raised.

    Parameters
    ----------
    organism : str
        KEGG organism code, e.g. ``"hsa"``, ``"mmu"``, ``"dre"``.
    output_dir : str
        Directory where the output files will be written.  Created if it does
        not exist.
    overwrite : bool
        If ``True``, existing files are overwritten.  If ``False`` (default),
        an error is raised when a target file already exists.
    filter_path : list of str or None
        Regular-expression patterns for pathways to remove
        (passed to ``buildGraphFromKEGGREST``).
    matrices : list of str or None
        Which matrices to compute.  Defaults to
        ``["hypergeom", "diffusion", "pagerank"]``
        (passed to ``buildDataFromGraph``).
    normality : list of str or None
        Methods for which row-sum vectors should be stored.  Defaults to
        ``["diffusion", "pagerank"]`` (passed to ``buildDataFromGraph``).
    damping_factor : float
        PageRank damping factor (default 0.85).
    niter : int
        Number of random trials for the null connected-component-size
        distribution (default 100, must be 10..1000).

    Returns
    -------
    dict
        Dictionary with keys ``"graph_path"``, ``"h5_path"``, ``"organism"``.

    Raises
    ------
    FileExistsError
        If a target file already exists and *overwrite* is ``False``.
    """
    # ------------------------------------------------------------------
    # Resolve paths
    # ------------------------------------------------------------------
    graph_filename = f"{organism}_graph.pkl"
    h5_filename = f"{organism}.h5"

    os.makedirs(output_dir, exist_ok=True)
    graph_path = os.path.join(output_dir, graph_filename)
    h5_path = os.path.join(output_dir, h5_filename)

    # ------------------------------------------------------------------
    # Check existing files
    # ------------------------------------------------------------------
    existing = []
    if os.path.exists(graph_path):
        existing.append(graph_path)
    if os.path.exists(h5_path):
        existing.append(h5_path)

    if existing and not overwrite:
        raise FileExistsError(
            f"Target file(s) already exist: {existing}. "
            f"Set overwrite=True to replace them."
        )

    # ------------------------------------------------------------------
    # 1. Build graph from KEGG REST
    # ------------------------------------------------------------------
    print(f"Building KEGG graph for organism: {organism}")
    g = buildGraphFromKEGGREST(organism=organism, filter_path=filter_path)

    # Save graph
    print(f"Saving graph to: {graph_path}")
    with open(graph_path, "wb") as f:
        pickle.dump(g, f)

    # ------------------------------------------------------------------
    # 2. Build HDF5 database from graph
    # ------------------------------------------------------------------
    print(f"Building HDF5 database: {h5_path}")
    buildDataFromGraph(
        keggdata_graph=g,
        h5_path=h5_path,
        matrices=matrices,
        normality=normality,
        damping_factor=damping_factor,
        niter=niter,
    )

    print("Dataset build complete.")
    return {
        "organism": organism,
        "graph_path": graph_path,
        "h5_path": h5_path,
    }
