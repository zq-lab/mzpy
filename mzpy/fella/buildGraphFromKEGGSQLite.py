"""
Python port of FELLA::buildGraphFromKEGGREST — SQLite backend

This module rebuilds the curated KEGG graph using a local SQLite database
instead of the KEGG REST API.  It expects a single table named ``item``
with at least the columns ``kid``, ``txt`` and ``category``.

The ``txt`` column should contain the raw KEGG flat-file text exactly as
returned by the KEGG API (e.g. ``get/C00001``).  The ``category`` column
classifies the entry (``pathway``, ``module``, ``enzyme``, ``reaction``,
``compound``, ``ko``, ``gene`` …).

The function signature and the returned ``igraph.Graph`` are kept identical
to ``buildGraphFromKEGGREST`` so that the two modules are drop-in replacements
for each other.

Dependencies
------------
- sqlite3  (stdlib)
- pandas
- igraph
"""

import os
import re
import sqlite3
from collections import defaultdict
from typing import Dict, List, Optional, Tuple

import pandas as pd
import igraph as ig

# ---------------------------------------------------------------------------
# KEGG flat-file parser
# ---------------------------------------------------------------------------

_KEGG_SECTION_LEN = 12


def _parse_kegg_sections(txt: str) -> Dict[str, List[str]]:
    """
    Parse a KEGG flat-file into a dict ``{section_name: [lines]}``.
    Continuation lines (12 leading spaces) are merged into the current
    section.  The terminator ``///`` stops parsing.
    """
    sections: Dict[str, List[str]] = defaultdict(list)
    current: Optional[str] = None
    for raw_line in txt.splitlines():
        line = raw_line.rstrip("\n")
        if line.startswith("///"):
            break
        if len(line) <= _KEGG_SECTION_LEN:
            continue
        field = line[:_KEGG_SECTION_LEN].rstrip()
        content = line[_KEGG_SECTION_LEN:].rstrip()
        if field:
            current = field
            sections[current].append(content)
        elif current is not None:
            sections[current].append(content)
    return dict(sections)


def _first_tokens(lines: List[str]) -> List[str]:
    """Return the first whitespace-separated token of every line."""
    out = []
    for line in lines:
        parts = line.split()
        if parts:
            out.append(parts[0])
    return out


def _extract_ecs_from_orthology(lines: List[str]) -> List[str]:
    """
    Extract EC numbers from ORTHOLOGY lines such as
    ``K00001  alcohol dehydrogenase [EC:1.1.1.1]``.
    """
    out = []
    for line in lines:
        # Find all [EC:...] blocks
        for m in re.finditer(r"\[EC:([\d.\-]+(?:\s+[\d.\-]+)*)\]", line):
            for ec in m.group(1).split():
                if "-" not in ec:
                    out.append(ec)
    return out


def _extract_ko_from_gene_lines(lines: List[str]) -> Dict[str, str]:
    """
    Parse GENE lines from an organism-specific pathway entry.
    ``10327  (K00001)  alcohol dehydrogenase`` → ``{gene_id: ko_id}``.
    The gene ID is returned *without* the organism prefix.
    """
    mapping = {}
    for line in lines:
        m = re.search(r"\(([Kk]\d+)\)", line)
        if m:
            parts = line.split()
            if parts:
                mapping[parts[0]] = m.group(1).upper()
    return mapping


def _extract_name(sections: Dict[str, List[str]]) -> str:
    """Best-effort extraction of the entry description / name."""
    for key in ("NAME", "DEFINITION", "TITLE"):
        if key in sections and sections[key]:
            return " ".join(sections[key])
    return ""


# ---------------------------------------------------------------------------
# ID sanitiser (mirrors FELLA:::sanitise)
# ---------------------------------------------------------------------------

def _sanitise(x: str, category: str, organism: str) -> Optional[str]:
    if category == "pathway":
        return re.sub(r"^path:(.+)(\d{5})$", r"\1\2", x)
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
    if category == "gene":
        return re.sub(rf"^{organism}:(.+)$", r"\1", x)
    if category == "ko":
        return re.sub(r"^ko:(K\d+)$", r"\1", x)
    return x


def _sanitise_ids(ids: List[str], category: str, organism: str) -> List[Optional[str]]:
    return [_sanitise(i, category, organism) for i in ids]


# ---------------------------------------------------------------------------
# SQLite equivalents of KEGG REST operations
# ---------------------------------------------------------------------------

def _sqlite_list(
    cursor: sqlite3.Cursor,
    category: str,
    organism: Optional[str] = None,
) -> Dict[str, str]:
    """
    Equivalent of ``keggList``.
    Returns ``{entry_id: definition}``.
    """
    if category == "pathway" and organism is not None:
        cursor.execute(
            """
            SELECT kid, txt FROM item
            WHERE category = 'pathway'
              AND (kid LIKE ? OR kid LIKE ?)
            """,
            (f"{organism}%", f"path:{organism}%"),
        )
    else:
        cursor.execute(
            "SELECT kid, txt FROM item WHERE category = ?",
            (category,),
        )
    out = {}
    for kid, txt in cursor.fetchall():
        sk = _sanitise(kid, category, organism or "")
        if sk is None:
            continue
        sections = _parse_kegg_sections(txt)
        out[sk] = _extract_name(sections)
    return out


def _sqlite_link(
    cursor: sqlite3.Cursor,
    target: str,
    source: str,
    organism: Optional[str] = None,
) -> Dict[str, List[str]]:
    """
    Equivalent of ``keggLink(target, source)``.
    Returns ``{source_id: [target_id, ...]}``.

    The implementation parses the raw KEGG flat-files stored in SQLite.
    """
    out: Dict[str, List[str]] = defaultdict(list)
    org = organism or ""

    # ------------------------------------------------------------------
    # Helper closures to keep the code readable
    # ------------------------------------------------------------------
    def _rows(cat: str, like_pat: Optional[Tuple[str, ...]] = None):
        if like_pat is not None:
            cursor.execute(
                """
                SELECT kid, txt FROM item
                WHERE category = ?
                  AND (kid LIKE ? OR kid LIKE ?)
                """,
                (cat, *like_pat),
            )
        else:
            cursor.execute(
                "SELECT kid, txt FROM item WHERE category = ?",
                (cat,),
            )
        return cursor.fetchall()

    def _sanitise_list(ids: List[str], cat: str) -> List[str]:
        return [s for s in _sanitise_ids(ids, cat, org) if s is not None]

    # ------------------------------------------------------------------
    # 1. pathway  as  source
    # ------------------------------------------------------------------
    if source == "pathway":
        like = (f"{org}%", f"path:{org}%") if organism else None
        for kid, txt in _rows("pathway", like):
            s_id = _sanitise(kid, "pathway", org)
            if s_id is None:
                continue
            sections = _parse_kegg_sections(txt)
            if target == "module":
                # Module entries contain PATHWAY section (reverse lookup)
                pass  # handled below when module is source
            elif target == "reaction":
                for t_id in _sanitise_list(_first_tokens(sections.get("REACTION", [])), "reaction"):
                    out[s_id].append(t_id)
            elif target == "compound":
                for t_id in _sanitise_list(_first_tokens(sections.get("COMPOUND", [])), "compound"):
                    out[s_id].append(t_id)
            elif target == "enzyme":
                # Direct enzyme links in pathways are not extracted here;
                # they are inferred via genes below.
                pass

    # ------------------------------------------------------------------
    # 2. module  as  source
    # ------------------------------------------------------------------
    if source == "module":
        for kid, txt in _rows("module"):
            s_id = _sanitise(kid, "module", org)
            if s_id is None:
                continue
            sections = _parse_kegg_sections(txt)
            if target == "pathway":
                for t_id in _sanitise_list(_first_tokens(sections.get("PATHWAY", [])), "pathway"):
                    out[s_id].append(t_id)
            elif target == "enzyme":
                for t_id in _sanitise_list(_extract_ecs_from_orthology(sections.get("ORTHOLOGY", [])), "enzyme"):
                    out[s_id].append(t_id)
            elif target == "reaction":
                for t_id in _sanitise_list(_first_tokens(sections.get("REACTION", [])), "reaction"):
                    out[s_id].append(t_id)
            elif target == "compound":
                for t_id in _sanitise_list(_first_tokens(sections.get("COMPOUND", [])), "compound"):
                    out[s_id].append(t_id)

    # ------------------------------------------------------------------
    # 3. enzyme  as  source
    # ------------------------------------------------------------------
    if source == "enzyme":
        for kid, txt in _rows("enzyme"):
            s_id = _sanitise(kid, "enzyme", org)
            if s_id is None:
                continue
            sections = _parse_kegg_sections(txt)
            if target == "reaction":
                for t_id in _sanitise_list(_first_tokens(sections.get("REACTION", [])), "reaction"):
                    out[s_id].append(t_id)
            elif target == "compound":
                for t_id in _sanitise_list(_first_tokens(sections.get("COMPOUND", [])), "compound"):
                    out[s_id].append(t_id)

    # ------------------------------------------------------------------
    # 4. reaction  as  source
    # ------------------------------------------------------------------
    if source == "reaction":
        for kid, txt in _rows("reaction"):
            s_id = _sanitise(kid, "reaction", org)
            if s_id is None:
                continue
            sections = _parse_kegg_sections(txt)
            if target == "compound":
                # Parse EQUATION to pull compound IDs (C#####)
                eq_ids = set()
                for line in sections.get("EQUATION", []):
                    eq_ids.update(re.findall(r"C\d{5}", line))
                for t_id in _sanitise_list(sorted(eq_ids), "compound"):
                    out[s_id].append(t_id)
            elif target == "enzyme":
                for t_id in _sanitise_list(_first_tokens(sections.get("ENZYME", [])), "enzyme"):
                    out[s_id].append(t_id)

    # ------------------------------------------------------------------
    # 5. compound  as  source
    # ------------------------------------------------------------------
    if source == "compound":
        for kid, txt in _rows("compound"):
            s_id = _sanitise(kid, "compound", org)
            if s_id is None:
                continue
            sections = _parse_kegg_sections(txt)
            if target == "reaction":
                for t_id in _sanitise_list(_first_tokens(sections.get("REACTION", [])), "reaction"):
                    out[s_id].append(t_id)
            elif target == "enzyme":
                for t_id in _sanitise_list(_first_tokens(sections.get("ENZYME", [])), "enzyme"):
                    out[s_id].append(t_id)
            elif target == "pathway":
                for t_id in _sanitise_list(_first_tokens(sections.get("PATHWAY", [])), "pathway"):
                    out[s_id].append(t_id)

    return dict(out)


def _sqlite_gene_mappings(
    cursor: sqlite3.Cursor,
    organism: str,
) -> Tuple[Dict[str, List[str]], Dict[str, List[str]], Dict[str, List[str]], Dict[str, List[str]]]:
    """
    Extract gene-centric mappings from the local SQLite.

    Returns
    -------
    m_path_gene : {pathway: [gene, ...]}
    m_mod_gene  : {module: [gene, ...]}
    m_gene_enzyme : {gene: [enzyme, ...]}
    kegg_gene2entrez : {gene: [entrez, ...]}

    Organism-specific pathways are parsed for their ``GENE`` section.
    KO entries are parsed for EC numbers in ``DEFINITION``.
    """
    m_path_gene: Dict[str, List[str]] = defaultdict(list)
    m_mod_gene: Dict[str, List[str]] = defaultdict(list)
    m_gene_enzyme: Dict[str, List[str]] = defaultdict(list)
    kegg_gene2entrez: Dict[str, List[str]] = defaultdict(list)

    # ------------------------------------------------------------------
    # 1. KO → EC cache (from ko entries)
    # ------------------------------------------------------------------
    cursor.execute("SELECT kid, txt FROM item WHERE category = 'ko'")
    ko2ec: Dict[str, List[str]] = defaultdict(list)
    for kid, txt in cursor.fetchall():
        ko_id = _sanitise(kid, "ko", organism)
        if ko_id is None:
            continue
        sections = _parse_kegg_sections(txt)
        for line in sections.get("DEFINITION", []):
            for m in re.finditer(r"\[EC:([\d.\-]+(?:\s+[\d.\-]+)*)\]", line):
                for ec in m.group(1).split():
                    if "-" not in ec:
                        ko2ec[ko_id].append(ec)
        # Also check NAME for EC hints (rare)
        for line in sections.get("NAME", []):
            for m in re.finditer(r"\[EC:([\d.\-]+(?:\s+[\d.\-]+)*)\]", line):
                for ec in m.group(1).split():
                    if "-" not in ec and ec not in ko2ec[ko_id]:
                        ko2ec[ko_id].append(ec)

    # ------------------------------------------------------------------
    # 2. Organism pathways → genes (GENE section)
    # ------------------------------------------------------------------
    cursor.execute(
        """
        SELECT kid, txt FROM item
        WHERE category = 'pathway'
          AND (kid LIKE ? OR kid LIKE ?)
        """,
        (f"{organism}%", f"path:{organism}%"),
    )
    for kid, txt in cursor.fetchall():
        path_id = _sanitise(kid, "pathway", organism)
        if path_id is None:
            continue
        sections = _parse_kegg_sections(txt)
        gene2ko = _extract_ko_from_gene_lines(sections.get("GENE", []))
        for gene_raw, ko_id in gene2ko.items():
            gene_id = f"{organism}:{gene_raw}"
            m_path_gene[path_id].append(gene_id)
            for ec in ko2ec.get(ko_id, []):
                m_gene_enzyme[gene_id].append(ec)

    # ------------------------------------------------------------------
    # 3. Modules → genes (via ORTHOLOGY if organism genes not stored)
    # ------------------------------------------------------------------
    cursor.execute("SELECT kid, txt FROM item WHERE category = 'module'")
    for kid, txt in cursor.fetchall():
        mod_id = _sanitise(kid, "module", organism)
        if mod_id is None:
            continue
        sections = _parse_kegg_sections(txt)
        for line in sections.get("ORTHOLOGY", []):
            m = re.search(r"\(([Kk]\d+)\)", line)
            if m:
                ko_id = m.group(1).upper()
                # Module-gene link is indirect; we store KO as a pseudo-gene
                # because the original inference only needs gene→enzyme.
                pseudo_gene = f"{organism}:{ko_id}"
                m_mod_gene[mod_id].append(pseudo_gene)
                for ec in ko2ec.get(ko_id, []):
                    m_gene_enzyme[pseudo_gene].append(ec)

    # ------------------------------------------------------------------
    # 4. NCBI gene IDs from DBLINKS (optional)
    # ------------------------------------------------------------------
    cursor.execute(
        """
        SELECT kid, txt FROM item
        WHERE category = 'gene'
          AND (kid LIKE ? OR kid LIKE ?)
        """,
        (f"{organism}%", f"{organism}:%"),
    )
    for kid, txt in cursor.fetchall():
        gene_id = _sanitise(kid, "gene", organism)
        if gene_id is None:
            continue
        sections = _parse_kegg_sections(txt)
        for line in sections.get("DBLINKS", []):
            if "NCBI-GeneID" in line or "GeneID" in line:
                # e.g. "NCBI-GeneID: 10327" or "GeneID: 10327"
                m = re.search(r":\s*(\d+)", line)
                if m:
                    kegg_gene2entrez[gene_id].append(m.group(1))

    return dict(m_path_gene), dict(m_mod_gene), dict(m_gene_enzyme), dict(kegg_gene2entrez)


# ---------------------------------------------------------------------------
# Inference helper (identical to REST version)
# ---------------------------------------------------------------------------

def _infere_con2ec(
    ids: List[str],
    ent: str,
    ent2gene: Dict[str, List[str]],
    gene2enzyme: Dict[str, List[str]],
) -> pd.DataFrame:
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
    df.attrs["from"] = "enzyme"
    df.attrs["to"] = ent
    return df


# ---------------------------------------------------------------------------
# Graph helpers (identical to REST version)
# ---------------------------------------------------------------------------

def _largestcc(graph: ig.Graph) -> ig.Graph:
    clusters = graph.clusters(mode="weak")
    return clusters.giant()


# ---------------------------------------------------------------------------
# Main function
# ---------------------------------------------------------------------------

def buildGraphFromKEGGSQLite(
    sqlite_path: str,
    organism: str = "hsa",
    filter_path: Optional[List[str]] = None,
) -> ig.Graph:
    """
    Build and return the curated KEGG graph from a local SQLite database.

    Parameters
    ----------
    sqlite_path : str
        Path to the SQLite file containing the ``item`` table.
    organism : str
        KEGG code for the organism of interest (default ``"hsa"``).
    filter_path : list of str or None
        Regular-expression patterns matching pathways that should be removed.

    Returns
    -------
    igraph.Graph
        Directed, curated KEGG graph with vertex attributes ``com``, ``NAME``,
        ``entrez`` and edge attribute ``weight``.
    """
    if not os.path.isfile(sqlite_path):
        raise FileNotFoundError(f"SQLite database not found: {sqlite_path}")

    categories = ["pathway", "module", "enzyme", "reaction", "compound"]
    conn = sqlite3.connect(sqlite_path)
    cursor = conn.cursor()

    # Verify table exists
    cursor.execute(
        "SELECT name FROM sqlite_master WHERE type='table' AND name='item'"
    )
    if not cursor.fetchone():
        conn.close()
        raise ValueError(
            f"Database {sqlite_path} does not contain the required 'item' table."
        )

    print(f"Building through local SQLite: {sqlite_path}")

    # ------------------------------------------------------------------
    # 1. List of id -> name for each category
    # ------------------------------------------------------------------
    print("Step 1/4: Loading KEGG entries from SQLite...")
    list_list: Dict[str, Dict[str, str]] = {}
    for category in categories:
        print(f"  -> category '{category}'...")
        if category == "pathway":
            raw = _sqlite_list(cursor, category, organism=organism)
        else:
            raw = _sqlite_list(cursor, category)
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
    # 2. Direct links (edges) between categories
    # ------------------------------------------------------------------
    print("Step 2/4: Extracting direct links from SQLite...")
    link_frames = []
    for j in range(len(categories)):          # source index
        for i in range(len(categories)):      # target index
            if i > j:
                tgt_db = categories[i]
                src_db = categories[j]
                print(f"  -> link/{tgt_db}/{src_db}...")
                raw_link = _sqlite_link(cursor, tgt_db, src_db, organism=organism)
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
    # 3. Gene-centric mappings (pathway/module -> enzyme via genes)
    # ------------------------------------------------------------------
    print("Step 3/4: Extracting gene-centric mappings from SQLite...")
    m_path_gene, m_mod_gene, m_gene_enzyme, kegg_gene2entrez = _sqlite_gene_mappings(
        cursor, organism
    )

    # enzyme -> organism genes (reverse of m_gene_enzyme)
    m_enzyme_gene: Dict[str, List[str]] = defaultdict(list)
    for gene, enzymes in m_gene_enzyme.items():
        for e in enzymes:
            m_enzyme_gene[e].append(gene)
    m_enzyme_gene = dict(m_enzyme_gene)

    # Map kegg enzymes to entrez (unique & sorted)
    m_enzyme_gene_entrez = {}
    for enz, genes in m_enzyme_gene.items():
        ncbis = set()
        for g in genes:
            for n in kegg_gene2entrez.get(g, []):
                ncbis.add(n)
        m_enzyme_gene_entrez[enz] = sorted(list(ncbis))

    print("  Done.")

    # ------------------------------------------------------------------
    # 4. Inferred connections (pathway/module -> enzyme via genes)
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
    # 5. Assemble edge list (identical logic to REST version)
    # ------------------------------------------------------------------
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
    df_edges = df_edges.dropna(subset=["from", "to"]).reset_index(drop=True)

    print("Done.")

    # ------------------------------------------------------------------
    # 6. Build raw graph (identical logic to REST version)
    # ------------------------------------------------------------------
    print("Building graph...")
    edge_tuples = list(zip(df_edges["from"], df_edges["to"]))
    if not edge_tuples:
        conn.close()
        raise RuntimeError("No edges were extracted from SQLite. Cannot build graph.")

    g_raw = ig.Graph.TupleList(edge_tuples, directed=True)
    g_raw = g_raw.simplify()

    coms = []
    for v in g_raw.vs["name"]:
        cat = map_category.get(v)
        coms.append(categories.index(cat) + 1 if cat is not None else None)
    g_raw.vs["com"] = coms

    missing_com = [i for i, c in enumerate(g_raw.vs["com"]) if c is None]
    if missing_com:
        g_raw.delete_vertices(missing_com)

    infere_from = set(df_infere["from"].dropna().unique()) if not df_infere.empty else set()
    bad_enzymes = [
        i for i, (c, n) in enumerate(zip(g_raw.vs["com"], g_raw.vs["name"]))
        if c == 3 and n not in infere_from
    ]
    if bad_enzymes:
        g_raw.delete_vertices(bad_enzymes)

    org_modules = set(m_mod_gene.keys())
    bad_modules = [
        i for i, (c, n) in enumerate(zip(g_raw.vs["com"], g_raw.vs["name"]))
        if c == 2 and n not in org_modules
    ]
    if bad_modules:
        g_raw.delete_vertices(bad_modules)

    perm = sorted(
        range(g_raw.vcount()),
        key=lambda idx: (g_raw.vs[idx]["com"], g_raw.vs[idx]["name"]),
    )
    rank = [0] * g_raw.vcount()
    for new_pos, old_idx in enumerate(perm):
        rank[old_idx] = new_pos
    g_raw = g_raw.permute_vertices(rank)

    weights = []
    for e in g_raw.es:
        src_com = g_raw.vs[e.source]["com"]
        tgt_com = g_raw.vs[e.target]["com"]
        weights.append(abs(src_com - tgt_com))
    g_raw.es["weight"] = weights

    w3_sources = {e.source for e in g_raw.es if e["weight"] == 3}
    bad_reactions = [
        i for i, c in enumerate(g_raw.vs["com"])
        if c == 4 and i not in w3_sources
    ]
    if bad_reactions:
        g_raw.delete_vertices(bad_reactions)

    w1_sources = {e.source for e in g_raw.es if e["weight"] == 1}
    bad_compounds = [
        i for i, c in enumerate(g_raw.vs["com"])
        if c == 5 and i not in w1_sources
    ]
    if bad_compounds:
        g_raw.delete_vertices(bad_compounds)

    if filter_path is not None:
        names_path = [v["name"] for v in g_raw.vs if v["com"] == 1]
        names_out = set()
        for pat in filter_path:
            names_out.update(n for n in names_path if re.search(pat, n))
        if names_out:
            print(f"Filtering {len(names_out)} pathways.")
            g_raw.delete_vertices([v.index for v in g_raw.vs if v["name"] in names_out])

    g_raw = _largestcc(g_raw)
    print("Done.")

    # ------------------------------------------------------------------
    # 7. Prune redundant transitive edges (identical logic)
    # ------------------------------------------------------------------
    print("Pruning graph...")
    edges_by_weight = defaultdict(list)
    for e in g_raw.es:
        edges_by_weight[e["weight"]].append(e.index)

    sorted_weights = sorted(edges_by_weight.keys())
    if not sorted_weights:
        conn.close()
        raise RuntimeError("Graph has no edges after initial filtering.")

    first_w = sorted_weights[0]
    print(f"Current weight: {first_w} out of 4...")
    g_curated = g_raw.subgraph_edges(edges_by_weight[first_w], delete_vertices=False)

    for w in sorted_weights[1:]:
        print(f"Current weight: {w} out of 4...")
        dist_matrix = g_curated.distances(mode="out")
        list_edges = edges_by_weight[w]
        list_ends = [(g_raw.es[eid].source, g_raw.es[eid].target) for eid in list_edges]

        new_edges = []
        new_weights = []
        for eid, (src, tgt) in zip(list_edges, list_ends):
            d = dist_matrix[src][tgt]
            if d > w:
                new_edges.append((src, tgt))
                new_weights.append(g_raw.es[eid]["weight"])

        if new_edges:
            g_curated.add_edges(new_edges, attributes={"weight": new_weights})

    for e in g_curated.es:
        e["weight"] = 1.0 / e["weight"]

    # ------------------------------------------------------------------
    # 8. Attach vertex metadata
    # ------------------------------------------------------------------
    tmp = {}
    for category in categories:
        tmp.update(list_list[category])

    names_attr = []
    entrez_attr = []
    for v in g_curated.vs:
        name = v["name"]
        desc = tmp.get(name, "")
        names_attr.append(desc.split("; ") if "; " in desc else [desc])
        entrez_attr.append(m_enzyme_gene_entrez.get(name, []))

    g_curated.vs["NAME"] = names_attr
    g_curated.vs["entrez"] = entrez_attr

    g_curated["organism"] = organism
    g_curated["comment"] = (
        f"SQLite backend; organism={organism}; "
        f"built from {sqlite_path}"
    )

    conn.close()
    print("Done.")
    return g_curated
