"""
Python port of FELLA::buildGraphFromKEGGREST

This module rebuilds the curated KEGG graph using the KEGG REST API.
It tries to stay as close as possible to the original R implementation
in terms of argument names, logic and edge curation steps.

Dependencies
------------
- requests
- igraph  (pip install igraph)
- pandas  (for intermediate data frames)
"""

import hashlib
import os
import re
import time
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
import pandas as pd
import igraph as ig

# ---------------------------------------------------------------------------
# KEGG REST helpers
# ---------------------------------------------------------------------------

_KEGG_BASE = "http://rest.kegg.jp"
_CACHE_DIR = Path.home() / ".cache" / "fella_py"
_CACHE_TTL_SECONDS = 7 * 24 * 3600  # 7 days


def _cache_path(url: str) -> Path:
    """Return local cache file path for a URL."""
    h = hashlib.md5(url.encode("utf-8")).hexdigest()
    return _CACHE_DIR / f"{h}.txt"


def _kegg_request(
    endpoint: str,
    max_retries: int = 5,
    backoff_factor: float = 2.0,
    timeout: float = 300.0,
    use_cache: bool = True,
) -> str:
    """
    Execute a GET request against the KEGG REST API with caching,
    session reuse, and automatic retry on transient failures.
    """
    url = f"{_KEGG_BASE}/{endpoint}"
    cache_file = _cache_path(url)

    # 1. Try cache
    if use_cache and cache_file.exists():
        if time.time() - cache_file.stat().st_mtime < _CACHE_TTL_SECONDS:
            with open(cache_file, "r", encoding="utf-8") as fh:
                return fh.read()

    # 2. Session with retry strategy
    session = requests.Session()
    retry_strategy = Retry(
        total=max_retries,
        backoff_factor=backoff_factor,
        status_forcelist=[429, 500, 502, 503, 504],
        allowed_methods=["GET", "POST"],
        raise_on_status=False,
    )
    adapter = HTTPAdapter(max_retries=retry_strategy)
    session.mount("http://", adapter)
    session.mount("https://", adapter)

    try:
        print(f"  [KEGG] GET {url}  (timeout={timeout}s, retries={max_retries})")
        resp = session.get(url, timeout=timeout)
        resp.raise_for_status()
        text = resp.text

        # 3. Save cache
        if use_cache:
            _CACHE_DIR.mkdir(parents=True, exist_ok=True)
            with open(cache_file, "w", encoding="utf-8") as fh:
                fh.write(text)
        return text
    except requests.RequestException as exc:
        raise RuntimeError(
            f"KEGG request failed after {max_retries} attempts: {url}"
        ) from exc
    finally:
        session.close()


def _kegg_info(organism: str) -> str:
    """Wrap keggInfo."""
    return _kegg_request(f"info/{organism}")


def _kegg_list(database: str, organism: Optional[str] = None) -> Dict[str, str]:
    """
    Wrap keggList.
    Returns {entry_id: definition}.
    Only pathways are organism-specific in the original code.
    """
    if organism is not None and database == "pathway":
        text = _kegg_request(f"list/{database}/{organism}")
    else:
        text = _kegg_request(f"list/{database}")
    out = {}
    for line in text.strip().splitlines():
        if "\t" in line:
            entry, definition = line.split("\t", 1)
            out[entry] = definition
    return out


def _kegg_link(target: str, source: str) -> Dict[str, List[str]]:
    """
    Wrap keggLink(target, source).
    Returns {source_id: [target_id, ...]}.
    """
    text = _kegg_request(f"link/{target}/{source}")
    out = defaultdict(list)
    for line in text.strip().splitlines():
        if "\t" in line:
            src, tgt = line.split("\t", 1)
            out[src].append(tgt)
    return out


def _kegg_conv(target_db: str, source_db: str) -> Dict[str, List[str]]:
    """
    Wrap keggConv(target_db, source_db).
    Returns {source_id: [target_id, ...]}.
    """
    text = _kegg_request(f"conv/{target_db}/{source_db}")
    out = defaultdict(list)
    for line in text.strip().splitlines():
        if "\t" in line:
            src, tgt = line.split("\t", 1)
            out[src].append(tgt)
    return out


# ---------------------------------------------------------------------------
# ID sanitiser
# ---------------------------------------------------------------------------

def _sanitise(x: str, category: str, organism: str) -> Optional[str]:
    """
    Sanitise a single KEGG identifier.
    Mirrors FELLA:::sanitise.
    """
    if category == "pathway":
        return re.sub(r"^path:(.+)(\d{5})$", f"{organism}\\2", x)
    if category == "module":
        return re.sub(r"^md:(.*)(M\d{5})$", r"\2", x)
    if category == "enzyme":
        if "-" in x:
            return None
        return re.sub(r"^ec:(\d+\.\d+\.\d+\.\d+)$", r"\1", x)
    if category == "ncbi":
        return re.sub(r"^ncbi[^:]+:(.*\d+)$", r"\1", x)
    if category == "reaction":
        return re.sub(r"^rn:(R\d{5})$", r"\1", x)
    if category == "compound":
        return re.sub(r"^cpd:(C\d{5})$", r"\1", x)
    return x


def _sanitise_ids(ids: List[str], category: str, organism: str) -> List[Optional[str]]:
    """Vectorised sanitise."""
    return [_sanitise(i, category, organism) for i in ids]


# ---------------------------------------------------------------------------
# Inference helper
# ---------------------------------------------------------------------------

def _infere_con2ec(
    ids: List[str],
    ent: str,
    ent2gene: Dict[str, List[str]],
    gene2enzyme: Dict[str, List[str]],
) -> pd.DataFrame:
    """
    Infer connections to EC families by passing through genes.
    Mirrors FELLA:::infere.con2ec.
    """
    rows = []
    for x in ids:
        genes = ent2gene.get(x, [])
        enzymes = set()
        for g in genes:
            for e in gene2enzyme.get(g, []):
                enzymes.add(e)
        for e in enzymes:
            rows.append({"from": e, "to": x})
    df = pd.DataFrame(rows)
    # store metadata (not used by pandas itself, but kept for parity)
    df.attrs["from"] = "enzyme"
    df.attrs["to"] = ent
    return df


# ---------------------------------------------------------------------------
# Graph helpers
# ---------------------------------------------------------------------------

def _largestcc(graph: ig.Graph) -> ig.Graph:
    """Extract the largest (weakly) connected component."""
    clusters = graph.clusters(mode="weak")
    giant = clusters.giant()
    return giant


# ---------------------------------------------------------------------------
# Main function
# ---------------------------------------------------------------------------

def buildGraphFromKEGGREST(
    organism: str = "hsa",
    filter_path: Optional[List[str]] = None,
) -> ig.Graph:
    """
    Build and return the curated KEGG graph via the KEGG REST API.

    Parameters
    ----------
    organism : str
        KEGG code for the organism of interest (default "hsa").
    filter_path : list of str or None
        Regular-expression patterns matching pathways that should be removed.
        E.g. ["01100"] filters the overview metabolic pathway.

    Returns
    -------
    igraph.Graph
        Directed, curated KEGG graph with vertex attributes 'com', 'NAME',
        'entrez' and edge attribute 'weight'.
    """
    categories = ["pathway", "module", "enzyme", "reaction", "compound"]

    # ------------------------------------------------------------------
    # 1. Organism info & gene annotation category
    # ------------------------------------------------------------------
    print("Building through KEGGREST...")
    info_org = _kegg_info(organism)

    # find available ncbi annotations
    info_geneannot = [line.strip().replace(" ", "") for line in info_org.splitlines()
                      if re.search(r"ncbi-[a-z]+", line)]
    cat_geneannot = next(
        (g for g in ("ncbi-geneid", "ncbi-proteinid") if g in info_geneannot),
        None,
    )
    if cat_geneannot is None:
        raise ValueError(
            f"Organism {organism} does not appear to have either "
            f"ncbi-geneid or ncbi-proteinid."
        )
    print(
        f"Available gene annotations: {', '.join(info_geneannot)}. "
        f"Using {cat_geneannot}"
    )

    # ------------------------------------------------------------------
    # 2. List of id -> name for each category
    # ------------------------------------------------------------------
    print("Step 1/4: Fetching KEGG entries for 5 categories...")
    list_list: Dict[str, Dict[str, str]] = {}
    for category in categories:
        print(f"  -> list/{category}...")
        if category == "pathway":
            raw = _kegg_list(category, organism=organism)
        else:
            raw = _kegg_list(category)
        # sanitise keys
        cleaned = {}
        for k, v in raw.items():
            sk = _sanitise(k, category, organism)
            if sk is not None:
                cleaned[sk] = v
        list_list[category] = cleaned
    print("  Done.")

    # map id -> category
    map_category = {}
    for category in categories:
        for entry_id in list_list[category]:
            map_category[entry_id] = category

    # ------------------------------------------------------------------
    # 3. Direct links (edges) between categories
    # ------------------------------------------------------------------
    print("Step 2/4: Fetching direct links between categories...")
    # Original R code uses lower.tri of a 5x5 matrix built from
    # expand.grid(categories, categories).  This selects target>source
    # combinations (by category index).
    link_frames = []
    for j in range(len(categories)):          # source index (Var2)
        for i in range(len(categories)):      # target index (Var1)
            if i > j:
                tgt_db = categories[i]
                src_db = categories[j]
                print(f"  -> link/{tgt_db}/{src_db}...")
                raw_link = _kegg_link(tgt_db, src_db)
                if not raw_link:
                    continue
                sources = []
                targets = []
                for src_id, tgt_ids in raw_link.items():
                    for tgt_id in tgt_ids:
                        sources.append(src_id)
                        targets.append(tgt_id)
                df = pd.DataFrame({
                    "from": _sanitise_ids(targets, tgt_db, organism),
                    "to": _sanitise_ids(sources, src_db, organism),
                })
                df.attrs["from"] = tgt_db
                df.attrs["to"] = src_db
                link_frames.append(df)
    print("  Done.")

    # ------------------------------------------------------------------
    # 4. Mappings through genes (for inferring EC connections)
    # ------------------------------------------------------------------
    print("Step 3/4: Fetching gene-centric mappings...")
    # pathway -> organism genes
    print("  -> link/pathway genes...")
    m_path_gene = _kegg_link(organism, "pathway")
    m_path_gene = {
        _sanitise(k, "pathway", organism): [_sanitise(v, "gene", organism) for v in vals]
        for k, vals in m_path_gene.items()
    }

    # module -> organism genes
    print("  -> link/module genes...")
    m_mod_gene = _kegg_link(organism, "module")
    m_mod_gene = {
        _sanitise(k, "module", organism): [_sanitise(v, "gene", organism) for v in vals]
        for k, vals in m_mod_gene.items()
    }

    # gene -> enzyme  (target=enzyme, source=organism)
    print("  -> link/enzyme genes...")
    m_gene_enzyme_raw = _kegg_link("enzyme", organism)
    m_gene_enzyme = defaultdict(list)
    for gene, enzymes in m_gene_enzyme_raw.items():
        sg = _sanitise(gene, "gene", organism)
        for e in enzymes:
            se = _sanitise(e, "enzyme", organism)
            if se is not None:
                m_gene_enzyme[sg].append(se)

    # kegg gene -> entrez (ncbi) ids
    print("  -> conv/gene IDs...")
    kegg_gene2entrez_raw = _kegg_conv(cat_geneannot, organism)
    kegg_gene2entrez = defaultdict(list)
    for gene, entrezs in kegg_gene2entrez_raw.items():
        sg = _sanitise(gene, "gene", organism)
        for e in entrezs:
            se = _sanitise(e, "ncbi", organism)
            if se is not None:
                kegg_gene2entrez[sg].append(se)

    # enzyme -> organism genes  (target=organism, source=enzyme)
    print("  -> link/enzyme reverse...")
    m_enzyme_gene_raw = _kegg_link(organism, "enzyme")
    m_enzyme_gene = defaultdict(list)
    for enz, genes in m_enzyme_gene_raw.items():
        se = _sanitise(enz, "enzyme", organism)
        for g in genes:
            sg = _sanitise(g, "gene", organism)
            if sg is not None and se is not None:
                m_enzyme_gene[se].append(sg)
    print("  Done.")

    # Map kegg enzymes to entrez (unique & sorted)
    m_enzyme_gene_entrez = {}
    for enz, genes in m_enzyme_gene.items():
        ncbis = set()
        for g in genes:
            for n in kegg_gene2entrez.get(g, []):
                ncbis.add(n)
        m_enzyme_gene_entrez[enz] = sorted(list(ncbis))

    # ------------------------------------------------------------------
    # 5. Inferred connections (pathway/module -> enzyme via genes)
    # ------------------------------------------------------------------
    con_infere = [
        _infere_con2ec(
            list(list_list["pathway"].keys()),
            "pathway",
            m_path_gene,
            m_gene_enzyme,
        ),
        _infere_con2ec(
            list(list_list["module"].keys()),
            "module",
            m_mod_gene,
            m_gene_enzyme,
        ),
    ]

    # ------------------------------------------------------------------
    # 6. Assemble edge list
    # ------------------------------------------------------------------
    # Direct connections: drop enzyme -> module/pathway because these are
    # inferred through genes instead.
    noinfere_frames = []
    for df in link_frames:
        a_from = df.attrs.get("from")
        a_to = df.attrs.get("to")
        if a_from == "enzyme" and a_to in ("module", "pathway"):
            continue
        noinfere_frames.append(df)

    df_noinfere = pd.concat(noinfere_frames, ignore_index=True) if noinfere_frames else pd.DataFrame(columns=["from", "to"])
    df_infere = pd.concat(con_infere, ignore_index=True) if con_infere else pd.DataFrame(columns=["from", "to"])

    df_edges = pd.concat([df_noinfere, df_infere], ignore_index=True)
    # drop rows with None ids (sanitise can produce NA for incomplete ECs)
    df_edges = df_edges.dropna(subset=["from", "to"]).reset_index(drop=True)

    print("Done.")

    # ------------------------------------------------------------------
    # 7. Build raw graph
    # ------------------------------------------------------------------
    print("Building graph...")
    edge_tuples = list(zip(df_edges["from"], df_edges["to"]))
    if not edge_tuples:
        raise RuntimeError("No edges were retrieved from KEGG. Cannot build graph.")

    g_raw = ig.Graph.TupleList(edge_tuples, directed=True)
    g_raw = g_raw.simplify()

    # assign community / category index (1-based, matching R)
    coms = []
    for v in g_raw.vs["name"]:
        cat = map_category.get(v)
        coms.append(categories.index(cat) + 1 if cat is not None else None)
    g_raw.vs["com"] = coms

    # delete vertices without a known category (obsolete / inexistent)
    missing_com = [i for i, c in enumerate(g_raw.vs["com"]) if c is None]
    if missing_com:
        g_raw.delete_vertices(missing_com)

    # delete enzymes not found in the inferred set (not in desired species)
    infere_from = set(df_infere["from"].dropna().unique()) if not df_infere.empty else set()
    bad_enzymes = [
        i for i, (c, n) in enumerate(zip(g_raw.vs["com"], g_raw.vs["name"]))
        if c == 3 and n not in infere_from
    ]
    if bad_enzymes:
        g_raw.delete_vertices(bad_enzymes)

    # delete modules that do not belong to the species
    org_modules = set(m_mod_gene.keys())
    bad_modules = [
        i for i, (c, n) in enumerate(zip(g_raw.vs["com"], g_raw.vs["name"]))
        if c == 2 and n not in org_modules
    ]
    if bad_modules:
        g_raw.delete_vertices(bad_modules)

    # order by category and id  (permute.vertices equivalent)
    perm = sorted(
        range(g_raw.vcount()),
        key=lambda idx: (g_raw.vs[idx]["com"], g_raw.vs[idx]["name"]),
    )
    rank = [0] * g_raw.vcount()
    for new_pos, old_idx in enumerate(perm):
        rank[old_idx] = new_pos
    g_raw = g_raw.permute_vertices(rank)

    # weight edges by absolute difference of category levels
    weights = []
    for e in g_raw.es:
        src_com = g_raw.vs[e.source]["com"]
        tgt_com = g_raw.vs[e.target]["com"]
        weights.append(abs(src_com - tgt_com))
    g_raw.es["weight"] = weights

    # keep only reactions that have at least one 3-weight edge (pathway link)
    w3_sources = {e.source for e in g_raw.es if e["weight"] == 3}
    bad_reactions = [
        i for i, c in enumerate(g_raw.vs["com"])
        if c == 4 and i not in w3_sources
    ]
    if bad_reactions:
        g_raw.delete_vertices(bad_reactions)

    # keep only compounds that have at least one 1-weight edge (reaction link)
    w1_sources = {e.source for e in g_raw.es if e["weight"] == 1}
    bad_compounds = [
        i for i, c in enumerate(g_raw.vs["com"])
        if c == 5 and i not in w1_sources
    ]
    if bad_compounds:
        g_raw.delete_vertices(bad_compounds)

    # filter pathways by regexp patterns
    if filter_path is not None:
        names_path = [v["name"] for v in g_raw.vs if v["com"] == 1]
        names_out = set()
        for pat in filter_path:
            names_out.update(n for n in names_path if re.search(pat, n))
        if names_out:
            print(f"Filtering {len(names_out)} pathways.")
            g_raw.delete_vertices([v.index for v in g_raw.vs if v["name"] in names_out])

    # keep largest connected component
    g_raw = _largestcc(g_raw)

    print("Done.")

    # ------------------------------------------------------------------
    # 8. Prune redundant transitive edges
    # ------------------------------------------------------------------
    print("Pruning graph...")
    # split edge ids by weight
    edges_by_weight = defaultdict(list)
    for e in g_raw.es:
        edges_by_weight[e["weight"]].append(e.index)

    # start with weight=1 edges
    sorted_weights = sorted(edges_by_weight.keys())
    if not sorted_weights:
        raise RuntimeError("Graph has no edges after initial filtering.")

    first_w = sorted_weights[0]
    print(f"Current weight: {first_w} out of 4...")
    g_curated = g_raw.subgraph_edges(edges_by_weight[first_w], delete_vertices=False)

    for w in sorted_weights[1:]:
        print(f"Current weight: {w} out of 4...")
        # all-pairs shortest path distances in the current curated graph (out-mode)
        dist_matrix = g_curated.distances(mode="out")

        list_edges = edges_by_weight[w]
        list_ends = [(g_raw.es[eid].source, g_raw.es[eid].target) for eid in list_edges]

        new_edges = []
        new_weights = []
        for eid, (src, tgt) in zip(list_edges, list_ends):
            # if no path exists, distance is inf -> keep the edge
            d = dist_matrix[src][tgt]
            if d > w:
                new_edges.append((src, tgt))
                new_weights.append(g_raw.es[eid]["weight"])

        if new_edges:
            g_curated.add_edges(new_edges, attributes={"weight": new_weights})

    # invert final edge weights
    for e in g_curated.es:
        e["weight"] = 1.0 / e["weight"]

    # ------------------------------------------------------------------
    # 9. Attach vertex metadata
    # ------------------------------------------------------------------
    # flattened name list: id -> description
    tmp = {}
    for category in categories:
        tmp.update(list_list[category])

    names_attr = []
    entrez_attr = []
    for v in g_curated.vs:
        name = v["name"]
        desc = tmp.get(name, "")
        names_attr.append(desc.split("; ") if "; " in desc else [desc])
        # entrez only meaningful for enzymes
        entrez_attr.append(m_enzyme_gene_entrez.get(name, []))

    g_curated.vs["NAME"] = names_attr
    g_curated.vs["entrez"] = entrez_attr

    # attach organism info as graph-level attributes
    g_curated["organism"] = organism
    # python-igraph does not have a native 'comment' slot; use a graph attribute
    g_curated["comment"] = info_org

    print("Done.")
    return g_curated
