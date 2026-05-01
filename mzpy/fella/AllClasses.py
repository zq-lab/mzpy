"""
Python mirrors of the S4 classes defined in FELLA::AllClasses.

These dataclasses replicate the slot structure of:
- D.keggdata, D.hypergeom, D.diffusion, D.pagerank
- FELLA.DATA
- U.userinput, U.hypergeom, U.diffusion, U.pagerank
- FELLA.USER
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Any

import numpy as np
import igraph as ig
from scipy import sparse


# ---------------------------------------------------------------------------
# Data classes (FELLA.DATA hierarchy)
# ---------------------------------------------------------------------------

@dataclass
class D_keggdata:
    graph: Optional[ig.Graph] = None
    id2name: Dict[str, List[str]] = field(default_factory=dict)
    pvalues_size: Optional[np.ndarray] = None
    id: Dict[str, Dict[str, int]] = field(default_factory=dict)
    status: str = "empty"


@dataclass
class D_hypergeom:
    matrix: Optional[sparse.csr_matrix] = None
    row_names: List[str] = field(default_factory=list)
    col_names: List[str] = field(default_factory=list)


@dataclass
class D_diffusion:
    matrix: Optional[np.ndarray] = None
    rowSums: Optional[Any] = None
    squaredRowSums: Optional[Any] = None


@dataclass
class D_pagerank:
    matrix: Optional[np.ndarray] = None
    rowSums: Optional[Any] = None
    squaredRowSums: Optional[Any] = None


@dataclass
class FELLA_DATA:
    keggdata: D_keggdata = field(default_factory=D_keggdata)
    hypergeom: D_hypergeom = field(default_factory=D_hypergeom)
    diffusion: D_diffusion = field(default_factory=D_diffusion)
    pagerank: D_pagerank = field(default_factory=D_pagerank)


# ---------------------------------------------------------------------------
# User classes (FELLA.USER hierarchy)
# ---------------------------------------------------------------------------

@dataclass
class U_userinput:
    metabolites: List[str] = field(default_factory=list)
    metabolitesbackground: List[str] = field(default_factory=list)
    excluded: List[str] = field(default_factory=list)


@dataclass
class U_hypergeom:
    valid: Optional[bool] = None
    pvalues: Optional[Any] = None
    pathhits: Optional[Any] = None
    pathbackground: Optional[Any] = None
    nbackground: Optional[int] = None
    ninput: Optional[int] = None


@dataclass
class U_diffusion:
    valid: Optional[bool] = None
    pscores: Optional[Any] = None
    approx: str = ""
    niter: Optional[int] = None


@dataclass
class U_pagerank:
    valid: Optional[bool] = None
    pscores: Optional[Any] = None
    approx: str = ""
    niter: Optional[int] = None


@dataclass
class FELLA_USER:
    userinput: U_userinput = field(default_factory=U_userinput)
    hypergeom: U_hypergeom = field(default_factory=U_hypergeom)
    diffusion: U_diffusion = field(default_factory=U_diffusion)
    pagerank: U_pagerank = field(default_factory=U_pagerank)
