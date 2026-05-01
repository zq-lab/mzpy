"""
Python port of FELLA::buildDataFromGraph

This module takes the KEGG graph produced by buildGraphFromKEGGREST and
computes the internal matrices / vectors required for enrichment analysis.
Results are stored in a single HDF5 file where the top-level group name is
the organism KEGG code (e.g. "hsa", "dre", "mmu").

Dependencies
------------
- igraph  (pip install igraph)
- numpy
- scipy   (sparse matrices + linear algebra)
- h5py
"""

import os
import pickle
from typing import List, Optional

import numpy as np
import igraph as ig
from scipy import linalg as la
from scipy import sparse
import h5py


# ---------------------------------------------------------------------------
# HDF5 helpers
# ---------------------------------------------------------------------------

def _save_graph_to_h5(grp, graph: ig.Graph):
    """Pickle an igraph object and store as uint8 array in HDF5."""
    b = pickle.dumps(graph)
    grp.create_dataset("graph", data=np.frombuffer(b, dtype=np.uint8))


def _save_csr_to_h5(grp, mat: sparse.csr_matrix, row_names: List[str], col_names: List[str]):
    """Store a CSR sparse boolean matrix into an HDF5 group."""
    grp.create_dataset("data", data=mat.data.astype(np.int8))
    grp.create_dataset("indices", data=mat.indices)
    grp.create_dataset("indptr", data=mat.indptr)
    grp.create_dataset("shape", data=np.array(mat.shape, dtype=np.int64))
    # store names as variable-length strings
    dt = h5py.string_dtype()
    grp.create_dataset("row_names", data=np.array(row_names, dtype=object), dtype=dt)
    grp.create_dataset("col_names", data=np.array(col_names, dtype=object), dtype=dt)


# ---------------------------------------------------------------------------
# Main function
# ---------------------------------------------------------------------------

def buildDataFromGraph(
    keggdata_graph: ig.Graph,
    h5_path: str,
    matrices: Optional[List[str]] = None,
    normality: Optional[List[str]] = None,
    damping_factor: float = 0.85,
    niter: int = 100,
) -> bool:
    """
    Build internal data files from a curated KEGG graph.

    Parameters
    ----------
    keggdata_graph : igraph.Graph
        The directed KEGG graph returned by buildGraphFromKEGGREST
        or buildGraphFromSQLite.
    h5_path : str
        Path to the output HDF5 file.  The file is always created anew
        (any existing file is overwritten).  The organism code is read
        from ``keggdata_graph["organism"]`` and stored as the root
        attribute ``organism``.
    matrices : list of str or None
        Which heavy matrices to store. Any of: "hypergeom", "diffusion",
        "pagerank".  Default is all three.
    normality : list of str or None
        For which methods the row-sum vectors should be stored. Any of:
        "diffusion", "pagerank".  Default is both.
    damping_factor : float
        Damping factor for PageRank (default 0.85).
    niter : int
        Number of random trials to estimate the null distribution of the
        largest connected-component size (default 100, must be 10..1000).

    Returns
    -------
    bool
        True if the database was written successfully.
            
    HDF5文件的内部结构
    dre_db.h5
    ├── @organism = "dre"        (root 属性，不再是 group)
    ├── graph                     (pickled igraph)
    ├── pvalues_size
    ├── hypergeom/
    │   ├── data, indices, indptr, shape, row_names, col_names
    ├── diffusion/
    │   ├── matrix, rowSums, squaredRowSums
    └── pagerank/
        ├── matrix, rowSums, squaredRowSums
    """
    # ------------------------------------------------------------------
    # Defaults & checks
    # ------------------------------------------------------------------
    if matrices is None:
        matrices = ["hypergeom", "diffusion", "pagerank"]
    if normality is None:
        normality = ["diffusion", "pagerank"]

    if not (0.0 < damping_factor < 1.0):
        raise ValueError("damping_factor must be between 0 and 1 (exclusive).")
    if not isinstance(h5_path, str) or not h5_path:
        raise ValueError("h5_path must be a non-empty string.")

    graph = keggdata_graph
    if graph is None:
        raise ValueError("keggdata_graph cannot be None.")
    if not isinstance(graph, ig.Graph):
        raise ValueError("keggdata_graph is not a valid igraph object.")
    if not graph.is_connected(mode="weak"):
        raise ValueError("keggdata_graph is not connected!")
    if "com" not in graph.vertex_attributes():
        raise ValueError(
            "keggdata_graph is not a valid graph: "
            "'com' attribute is missing."
        )

    valid_methods = {"hypergeom", "diffusion", "pagerank"}
    if matrices is not None and not set(matrices).issubset(valid_methods):
        raise ValueError(
            "matrices should contain one or more of: "
            "'hypergeom', 'diffusion', 'pagerank'."
        )
    if normality is not None and not set(normality).issubset(
        {"diffusion", "pagerank"}
    ):
        raise ValueError(
            "normality should contain one or more of: "
            "'diffusion', 'pagerank'."
        )

    if not (10 <= niter <= 1000):
        raise ValueError("niter must be an integer between 10 and 1000.")
    niter = int(round(niter))

    # ------------------------------------------------------------------
    # Resolve organism
    # ------------------------------------------------------------------
    try:
        organism = graph["organism"]
    except KeyError:
        organism = "unknown"

    # Ensure parent directory exists
    parent = os.path.dirname(os.path.abspath(h5_path))
    if parent:
        os.makedirs(parent, exist_ok=True)

    # ------------------------------------------------------------------
    # 1. Null distribution of largest-CC sizes under random sampling
    # ------------------------------------------------------------------
    max_subgraph_size = 250
    max_component_size = 250
    subgraph_sizes = np.arange(1, max_subgraph_size + 1)
    component_sizes = np.arange(1, max_component_size + 1)

    print(
        "Computing probabilities for random subgraphs... "
        "(this may take a while)"
    )

    n_nodes = graph.vcount()
    # array shape: (niter, subgraph_size, component_size)
    array_null = np.zeros(
        (niter, max_subgraph_size, max_component_size), dtype=bool
    )

    for trial in range(niter):
        # sample max_subgraph_size distinct nodes once
        select = np.random.choice(n_nodes, size=max_subgraph_size, replace=False)
        for idx, k_sub in enumerate(subgraph_sizes):
            vids = select[:k_sub].tolist()
            subg = graph.induced_subgraph(vids)
            clusters = subg.clusters(mode="weak")
            csize_max = max(clusters.sizes()) if clusters.sizes() else 0
            array_null[trial, idx, :] = component_sizes <= csize_max

    # empirical p-value-like quantity; transpose to (component, subgraph)
    keggdata_pvalues_size = (
        (array_null.sum(axis=0) + 1) / (niter + 1)
    ).T

    # ------------------------------------------------------------------
    # Write HDF5 file (always overwrite)
    # ------------------------------------------------------------------
    with h5py.File(h5_path, "w") as h5f:
        # Store organism as a root attribute for easy identification
        h5f.attrs["organism"] = organism

        # Save graph + p-values
        _save_graph_to_h5(h5f, graph)
        h5f.create_dataset("pvalues_size", data=keggdata_pvalues_size)
        print("Done.")

        # vertex indices for the two extreme levels (0-based in python-igraph)
        id_pathway = [v.index for v in graph.vs if v["com"] == 1]
        id_compound = [v.index for v in graph.vs if v["com"] == 5]

        # ------------------------------------------------------------------
        # 2. Hypergeometric matrix
        # ------------------------------------------------------------------
        if "hypergeom" in matrices:
            print("Computing hypergeom_matrix...")

            # In the original graph weights are 1/level_diff.
            # Reverse them back to the raw level differences for shortest-path.
            g = graph.copy()
            for e in g.es:
                e["weight"] = 1.0 / e["weight"]

            dists = np.array(
                g.shortest_paths(
                    source=id_compound, target=id_pathway, mode="out"
                )
            )
            hypergeom_matrix = sparse.csr_matrix(dists == 4, dtype=bool)
            hypergeom_row_names = [v["name"] for v in graph.vs if v["com"] == 5]
            hypergeom_col_names = [v["name"] for v in graph.vs if v["com"] == 1]

            hg_grp = h5f.create_group("hypergeom")
            _save_csr_to_h5(hg_grp, hypergeom_matrix, hypergeom_row_names, hypergeom_col_names)
            print("Done.")

        # ------------------------------------------------------------------
        # 3. Diffusion files
        # ------------------------------------------------------------------
        if "diffusion" in matrices or "diffusion" in normality:
            print(
                "Computing diffusion_matrix... "
                "(this may take a while and use some memory)"
            )

            g = graph.as_undirected()
            # Laplacian (dense numpy array)
            K = np.array(g.laplacian(normalized=False))

            # Boundary conditions: add 1 to diagonal entries of pathway nodes
            diag = np.diag(K).copy()
            diag[id_pathway] += 1.0
            np.fill_diagonal(K, diag)

            # Inverse (dense, memory-intensive)
            R = la.inv(K)
            # diffusion.matrix = all rows x compound columns
            diffusion_matrix = R[:, np.ix_(id_compound)[0]]

            if "diffusion" in matrices:
                dm_grp = h5f.create_group("diffusion")
                dm_grp.create_dataset("matrix", data=diffusion_matrix)
                print("Done.")
            else:
                # still need group for rowSums
                dm_grp = h5f.create_group("diffusion")

            if "diffusion" in normality:
                print("Computing diffusion_rowSums...")
                diffusion_rowSums = diffusion_matrix.sum(axis=1)
                diffusion_squaredRowSums = np.square(diffusion_matrix).sum(axis=1)
                dm_grp.create_dataset("rowSums", data=diffusion_rowSums)
                dm_grp.create_dataset("squaredRowSums", data=diffusion_squaredRowSums)
                print("Done.")

        # ------------------------------------------------------------------
        # 4. PageRank files
        # ------------------------------------------------------------------
        if "pagerank" in matrices or "pagerank" in normality:
            print(
                "Computing pagerank_matrix... "
                "(this may take a while and use some memory)"
            )

            g = graph.copy()
            d = damping_factor

            # Laplacian of the directed graph, transposed and negated
            K = -np.array(g.laplacian(normalized=False)).T
            np.fill_diagonal(K, 0.0)

            # Normalise columns to stochastic vectors
            n = K.shape[0]
            for j in range(n):
                col_sum = K[:, j].sum()
                if col_sum != 0:
                    K[:, j] /= col_sum
                else:
                    K[:, j] = 1.0 / n

            # Fundamental matrix: (I - d * K)^-1 * (1 - d)
            R = (1 - d) * la.inv(np.eye(n) - d * K)
            pagerank_matrix = R[:, np.ix_(id_compound)[0]]

            if "pagerank" in matrices:
                pr_grp = h5f.create_group("pagerank")
                pr_grp.create_dataset("matrix", data=pagerank_matrix)
                print("Done.")
            else:
                pr_grp = h5f.create_group("pagerank")

            if "pagerank" in normality:
                print("Computing pagerank_rowSums...")
                pagerank_rowSums = pagerank_matrix.sum(axis=1)
                pagerank_squaredRowSums = np.square(pagerank_matrix).sum(axis=1)
                pr_grp.create_dataset("rowSums", data=pagerank_rowSums)
                pr_grp.create_dataset("squaredRowSums", data=pagerank_squaredRowSums)
                print("Done.")

    return True
