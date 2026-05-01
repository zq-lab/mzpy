"""
Python port of FELLA::loadKEGGdata

Loads the internal database files (HDF5) produced by buildDataFromGraph
into a Python object tree that mirrors the original S4 FELLA.DATA hierarchy.

Dependencies
------------
- igraph  (pip install igraph)
- numpy
- scipy
- pandas  (for named vectors)
- h5py
"""

import os
import pickle
from typing import List, Optional

import numpy as np
import pandas as pd
import igraph as ig
from scipy import sparse
import h5py

from .AllClasses import (
    FELLA_DATA,
    D_keggdata,
    D_hypergeom,
    D_diffusion,
    D_pagerank,
)


# ---------------------------------------------------------------------------
# HDF5 helpers
# ---------------------------------------------------------------------------

def _load_graph_from_h5(h5f):
    """Load a pickled igraph object from an HDF5 dataset."""
    b = h5f["graph"][()].tobytes()
    return pickle.loads(b)


def _load_csr_from_h5(h5f, group_name: str):
    """Reconstruct a CSR sparse matrix from an HDF5 sub-group."""
    grp = h5f[group_name]
    data = grp["data"][()]
    indices = grp["indices"][()]
    indptr = grp["indptr"][()]
    shape = tuple(grp["shape"][()])
    row_names = [s.decode("utf-8") if isinstance(s, bytes) else s for s in grp["row_names"][()]]
    col_names = [s.decode("utf-8") if isinstance(s, bytes) else s for s in grp["col_names"][()]]
    mat = sparse.csr_matrix((data, indices, indptr), shape=shape, dtype=bool)
    return mat, row_names, col_names


# ---------------------------------------------------------------------------
# Main function
# ---------------------------------------------------------------------------

def loadKEGGdata(
    h5_path: str,
    load_matrix: Optional[List[str]] = None,
) -> FELLA_DATA:
    """
    Load a previously built KEGG knowledge database into a FELLA_DATA object.

    Parameters
    ----------
    h5_path : str
        Path to the HDF5 file produced by ``buildDataFromGraph``.
    load_matrix : list of str or None
        Which heavy matrices to load into memory. Can contain any of:
        ``"diffusion"``, ``"pagerank"``.  If None, no heavy matrices are
        loaded (row-sum vectors are still loaded when available).

    Returns
    -------
    FELLA_DATA
        Object containing the KEGG graph, p-value tables, and optionally the
        hypergeometric / diffusion / pagerank matrices.
    """
    # ------------------------------------------------------------------
    # Argument checks
    # ------------------------------------------------------------------
    if not isinstance(h5_path, str):
        raise ValueError("h5_path must be a string.")
    if not os.path.exists(h5_path):
        raise FileNotFoundError(f"HDF5 database file not found: {h5_path}")

    if load_matrix is not None:
        if not isinstance(load_matrix, (list, tuple)):
            raise ValueError("load_matrix must be a list/tuple or None.")
        valid = {"diffusion", "pagerank"}
        if not set(load_matrix).issubset(valid):
            raise ValueError(
                "load_matrix can only contain 'diffusion' and/or 'pagerank'."
            )

    # ------------------------------------------------------------------
    # Open HDF5 and read from root (flat structure)
    # ------------------------------------------------------------------
    with h5py.File(h5_path, "r") as h5f:
        organism = h5f.attrs.get("organism", "unknown")
        print(f"Loading database for organism: {organism}")

        # ------------------------------------------------------------------
        # 1. Load graph & p-values (required)
        # ------------------------------------------------------------------
        print("Loading KEGG graph data...")
        graph = _load_graph_from_h5(h5f)
        pvalues_size = h5f["pvalues_size"][()] if "pvalues_size" in h5f else None

        f_data = FELLA_DATA()
        f_data.keggdata.graph = graph
        f_data.keggdata.pvalues_size = pvalues_size

        # id2name: KEGG id -> list of names
        f_data.keggdata.id2name = {
            v["name"]: v["NAME"] for v in graph.vs if "NAME" in v.attributes()
        }

        # id: category -> {name: vertex_index}
        for i, cat in enumerate(
            ["pathway", "module", "enzyme", "reaction", "compound"], start=1
        ):
            f_data.keggdata.id[cat] = {
                v["name"]: v.index for v in graph.vs if v.attributes().get("com") == i
            }
        print("Done.")

        # ------------------------------------------------------------------
        # 2. Hypergeometric matrix
        # ------------------------------------------------------------------
        print("Loading hypergeom data...")
        print("Loading matrix...")
        if "hypergeom" in h5f:
            mat, rnames, cnames = _load_csr_from_h5(h5f, "hypergeom")
            f_data.hypergeom.matrix = mat
            f_data.hypergeom.row_names = rnames
            f_data.hypergeom.col_names = cnames
        else:
            print(
                "'hypergeom' group not present in HDF5. "
                "Hypergeometric test won't execute."
            )
        print("Done.")

        # ------------------------------------------------------------------
        # 3. Diffusion files
        # ------------------------------------------------------------------
        print("Loading diffusion data...")
        print("Loading matrix...")
        if load_matrix is not None and "diffusion" in load_matrix:
            if "diffusion" in h5f:
                dm_grp = h5f["diffusion"]
                if "matrix" in dm_grp:
                    f_data.diffusion.matrix = dm_grp["matrix"][()]
                else:
                    print(
                        "'diffusion/matrix' not present in HDF5. "
                        "Simulated permutations may execute slower for diffusion."
                    )
            else:
                print(
                    "'diffusion' group not present in HDF5. "
                    "Simulated permutations may execute slower for diffusion."
                )
        else:
            print(
                "'diffusion/matrix' not loaded. "
                "Simulated permutations may execute slower for diffusion."
            )
        print("Done.")

        print("Loading rowSums...")
        if "diffusion" in h5f:
            dm_grp = h5f["diffusion"]
            if "rowSums" in dm_grp and "squaredRowSums" in dm_grp:
                names = graph.vs["name"]
                f_data.diffusion.rowSums = pd.Series(
                    dm_grp["rowSums"][()], index=names
                )
                f_data.diffusion.squaredRowSums = pd.Series(
                    dm_grp["squaredRowSums"][()], index=names
                )
            else:
                print(
                    "'diffusion' rowSums not present in HDF5. "
                    "Z-scores won't be available for diffusion."
                )
        else:
            print(
                "'diffusion' group not present in HDF5. "
                "Z-scores won't be available for diffusion."
            )
        print("Done.")

        # ------------------------------------------------------------------
        # 4. PageRank files
        # ------------------------------------------------------------------
        print("Loading pagerank data...")
        print("Loading matrix...")
        if load_matrix is not None and "pagerank" in load_matrix:
            if "pagerank" in h5f:
                pr_grp = h5f["pagerank"]
                if "matrix" in pr_grp:
                    f_data.pagerank.matrix = pr_grp["matrix"][()]
                else:
                    print(
                        "'pagerank/matrix' not present in HDF5. "
                        "Simulated permutations may execute slower for pagerank."
                    )
            else:
                print(
                    "'pagerank' group not present in HDF5. "
                    "Simulated permutations may execute slower for pagerank."
                )
        else:
            print(
                "'pagerank/matrix' not loaded. "
                "Simulated permutations may execute slower for pagerank."
            )
        print("Done.")

        print("Loading rowSums...")
        if "pagerank" in h5f:
            pr_grp = h5f["pagerank"]
            if "rowSums" in pr_grp and "squaredRowSums" in pr_grp:
                names = graph.vs["name"]
                f_data.pagerank.rowSums = pd.Series(
                    pr_grp["rowSums"][()], index=names
                )
                f_data.pagerank.squaredRowSums = pd.Series(
                    pr_grp["squaredRowSums"][()], index=names
                )
            else:
                print(
                    "'pagerank' rowSums not present in HDF5. "
                    "Z-scores won't be available for pagerank."
                )
        else:
            print(
                "'pagerank' group not present in HDF5. "
                "Z-scores won't be available for pagerank."
            )
        print("Done.")

    f_data.keggdata.status = "loaded"
    print("Data successfully loaded.")
    return f_data
