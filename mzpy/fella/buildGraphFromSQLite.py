"""
buildGraphFromSQLite
====================

从本地 SQLite 数据库的 ``item`` 表中读取 KEGG 原始文本，
构建与 ``buildGraphFromKEGGREST`` 输出完全一致的 ``igraph.Graph``。

表结构要求
----------
``item`` 表至少包含以下列：

- ``kid``       : KEGG ID，如 ``C00001``、``R00001``、``dre00010``、``1.1.1.1``
- ``txt``       : 通过 ``rest.kegg.jp/get/{kid}`` 获取的原始 flat-file 文本
- ``update_on`` : 数据更新时间（本函数不使用）
- ``category``  : 数据类型，如 ``compound``、``reaction``、``pathway``、
                  ``module``、``enzyme``、``ko``、``gene``

调用示例
--------
>>> from fella_py.buildGraphFromSQLite import buildGraphFromSQLite
>>> g = buildGraphFromSQLite("dre", ["01100"], "/path/to/kegg.db")
"""

import os
import re
import sqlite3
from collections import defaultdict
from typing import Dict, List, Optional

import pandas as pd
import igraph as ig

# ---------------------------------------------------------------------------
# KEGG flat-file 解析器（参考 kegg.py 中 Node._parse 的思想）
# ---------------------------------------------------------------------------

def _parse_entry(txt: str) -> Dict[str, str]:
    """
    把 KEGG flat-file 文本拆分为 {节名: 节内容} 字典。
    节名全部转为小写，方便后续统一访问。
    """
    sections: Dict[str, str] = {}
    if not txt or not txt.strip():
        return sections

    # 找出所有节标题（行首的单词，后接至少一个空格）
    tags = re.findall(r"(\A\w+ )", txt) + re.findall(r"(\n\w+ )", txt)
    tags.append("///")

    for i in range(len(tags) - 1):
        pattern = re.compile(
            re.escape(tags[i]) + r"(.*?)" + re.escape(tags[i + 1]),
            re.DOTALL,
        )
        result = pattern.search(txt)
        key = tags[i].strip().lower()
        if result:
            sections[key] = result.group(1).strip()
        else:
            sections[key] = ""
    return sections


def _extract_name(txt: str) -> str:
    """从 KEGG 条目中提取名称 / 定义。"""
    sections = _parse_entry(txt)
    if "name" in sections and sections["name"]:
        # 取第一行，去掉末尾分号
        first = sections["name"].split("\n", 1)[0].strip()
        return first.rstrip(";").strip()
    if "definition" in sections and sections["definition"]:
        return sections["definition"].split("\n", 1)[0].strip()
    return ""


def _extract_linked_ids(txt: str, target_type: str) -> List[str]:
    """
    从单条 KEGG 文本中提取与 ``target_type`` 关联的 KEGG ID。

    target_type 支持: pathway, module, enzyme, reaction, compound, ko。
    """
    sections = _parse_entry(txt)

    def _find_ids(pattern: str, *keys: str) -> List[str]:
        ids: set = set()
        for key in keys:
            if key in sections and sections[key]:
                ids.update(re.findall(pattern, sections[key]))
        return list(ids)

    if target_type == "compound":
        # compound / equation / substrate / product 中均可能出现 C#####
        return _find_ids(r"C\d{5}", "compound", "equation", "substrate", "product")

    if target_type == "reaction":
        return _find_ids(r"R\d{5}", "reaction")

    if target_type == "enzyme":
        ids: set = set()
        # 1) 直接从 enzyme 节提取完整 EC 号
        for key in ("enzyme",):
            if key in sections and sections[key]:
                ids.update(re.findall(r"\d+\.\d+\.\d+\.\d+", sections[key]))
        # 2) 从 [EC:...] 块中提取（常见于 ORTHOLOGY / GENE 节）
        for key in ("orthology", "gene"):
            if key in sections and sections[key]:
                for m in re.finditer(
                    r"\[EC:([\d.\-]+(?:\s+[\d.\-]+)*)\]", sections[key]
                ):
                    for ec in m.group(1).split():
                        if "-" not in ec:
                            ids.add(ec)
        return list(ids)

    if target_type == "module":
        return _find_ids(r"M\d{5}", "module")

    if target_type == "pathway":
        return _find_ids(r"[a-z]{3,4}\d{5}", "pathway")

    if target_type == "ko":
        return _find_ids(r"K\d{5}", "orthology", "ko")

    return []


# ---------------------------------------------------------------------------
# ID 清洗（与 buildGraphFromKEGGREST 保持一致，仅修复 compound 的 bug）
# ---------------------------------------------------------------------------

def _sanitise(x: str, category: str, organism: str) -> Optional[str]:
    if not x:
        return None
    if category == "pathway":
        return re.sub(r"^path:(.+)(\d{5})$", rf"{organism}\2", x)
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
        # 原版 REST 代码此处写成了 \2（无 group 2），修正为 \1
        return re.sub(r"^cpd:(C\d{5})$", r"\1", x)
    if category == "gene":
        return re.sub(rf"^{organism}:(.+)$", r"\1", x)
    if category == "ko":
        return re.sub(r"^ko:(K\d+)$", r"\1", x)
    return x


def _sanitise_ids(ids: List[str], category: str, organism: str) -> List[Optional[str]]:
    return [_sanitise(i, category, organism) for i in ids]


# ---------------------------------------------------------------------------
# SQLite 版 KEGG REST 等价操作
# ---------------------------------------------------------------------------

def _sqlite_list(
    cursor: sqlite3.Cursor,
    category: str,
    organism: Optional[str] = None,
) -> Dict[str, str]:
    """等价于 keggList，返回 {entry_id: definition}。"""
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

    out: Dict[str, str] = {}
    for kid, txt in cursor.fetchall():
        sk = _sanitise(kid, category, organism or "")
        if sk is not None:
            out[sk] = _extract_name(txt)
    return out


def _sqlite_link(
    cursor: sqlite3.Cursor,
    target: str,
    source: str,
    organism: Optional[str] = None,
) -> Dict[str, List[str]]:
    """
    等价于 keggLink(target, source)。
    返回 {source_id: [target_id, ...]}。
    """
    out: Dict[str, List[str]] = defaultdict(list)
    org = organism or ""

    # 确定要查询的 category 与过滤条件
    if source == "pathway":
        cursor.execute(
            """
            SELECT kid, txt FROM item
            WHERE category = 'pathway'
              AND (kid LIKE ? OR kid LIKE ?)
            """,
            (f"{org}%", f"path:{org}%"),
        )
    elif source == "module":
        cursor.execute("SELECT kid, txt FROM item WHERE category = 'module'")
    elif source == "enzyme":
        cursor.execute("SELECT kid, txt FROM item WHERE category = 'enzyme'")
    elif source == "reaction":
        cursor.execute("SELECT kid, txt FROM item WHERE category = 'reaction'")
    elif source == "compound":
        cursor.execute("SELECT kid, txt FROM item WHERE category = 'compound'")
    else:
        return dict(out)

    for kid, txt in cursor.fetchall():
        src_id = _sanitise(kid, source, org)
        if src_id is None:
            continue
        for raw_tid in _extract_linked_ids(txt, target):
            tid = _sanitise(raw_tid, target, org)
            if tid is not None:
                out[src_id].append(tid)

    return dict(out)


def _sqlite_gene_mappings(
    cursor: sqlite3.Cursor,
    organism: str,
) -> tuple:
    """
    从本地 SQLite 提取基因层面的映射，用于推断 pathway/module → enzyme。

    Returns
    -------
    m_path_gene, m_mod_gene, m_gene_enzyme, kegg_gene2entrez
    """
    m_path_gene: Dict[str, List[str]] = defaultdict(list)
    m_mod_gene: Dict[str, List[str]] = defaultdict(list)
    m_gene_enzyme: Dict[str, List[str]] = defaultdict(list)
    kegg_gene2entrez: Dict[str, List[str]] = defaultdict(list)

    # ------------------------------------------------------------------
    # 1. 缓存 KO → EC（从 ko 条目解析 DEFINITION / NAME 中的 [EC:...]）
    # ------------------------------------------------------------------
    ko2ec: Dict[str, List[str]] = {}
    cursor.execute("SELECT kid, txt FROM item WHERE category = 'ko'")
    for kid, txt in cursor.fetchall():
        ko_id = _sanitise(kid, "ko", organism)
        if not ko_id:
            continue
        sections = _parse_entry(txt)
        ec_set: set = set()
        for key in ("definition", "name"):
            if key in sections and sections[key]:
                for m in re.finditer(
                    r"\[EC:([\d.\-]+(?:\s+[\d.\-]+)*)\]", sections[key]
                ):
                    for ec in m.group(1).split():
                        if "-" not in ec:
                            ec_set.add(ec)
        ko2ec[ko_id] = sorted(ec_set)

    # ------------------------------------------------------------------
    # 2. organism pathway → gene（解析 GENE 节中的 gene_id 与 KO）
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
        if not path_id:
            continue
        sections = _parse_entry(txt)
        if "gene" not in sections or not sections["gene"]:
            continue
        for line in sections["gene"].split("\n"):
            line = line.strip()
            m_gene = re.match(r"^(\d+)", line)
            m_ko = re.search(r"\(([Kk]\d+)\)", line)
            if m_gene and m_ko:
                gene_id = f"{organism}:{m_gene.group(1)}"
                ko_id = m_ko.group(1).upper()
                m_path_gene[path_id].append(gene_id)
                for ec in ko2ec.get(ko_id, []):
                    m_gene_enzyme[gene_id].append(ec)

    # ------------------------------------------------------------------
    # 3. module → KO（作为 pseudo-gene）→ EC
    # ------------------------------------------------------------------
    cursor.execute("SELECT kid, txt FROM item WHERE category = 'module'")
    for kid, txt in cursor.fetchall():
        mod_id = _sanitise(kid, "module", organism)
        if not mod_id:
            continue
        sections = _parse_entry(txt)
        if "orthology" not in sections or not sections["orthology"]:
            continue
        for line in sections["orthology"].split("\n"):
            line = line.strip()
            m_ko = re.search(r"\(([Kk]\d+)\)", line)
            if m_ko:
                ko_id = m_ko.group(1).upper()
                pseudo_gene = f"{organism}:{ko_id}"
                m_mod_gene[mod_id].append(pseudo_gene)
                for ec in ko2ec.get(ko_id, []):
                    m_gene_enzyme[pseudo_gene].append(ec)

    # ------------------------------------------------------------------
    # 4. gene → NCBI GeneID（从 DBLINKS 节提取，可选）
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
        if not gene_id:
            continue
        sections = _parse_entry(txt)
        if "dblinks" not in sections or not sections["dblinks"]:
            continue
        for line in sections["dblinks"].split("\n"):
            if "NCBI-GeneID" in line or "GeneID" in line:
                m = re.search(r":\s*(\d+)", line)
                if m:
                    kegg_gene2entrez[gene_id].append(m.group(1))

    return (
        dict(m_path_gene),
        dict(m_mod_gene),
        dict(m_gene_enzyme),
        dict(kegg_gene2entrez),
    )


# ---------------------------------------------------------------------------
# Inference helper（与 REST 版完全一致）
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
# Graph helpers（与 REST 版完全一致）
# ---------------------------------------------------------------------------

def _largestcc(graph: ig.Graph) -> ig.Graph:
    clusters = graph.clusters(mode="weak")
    return clusters.giant()


# ---------------------------------------------------------------------------
# 主函数
# ---------------------------------------------------------------------------

def buildGraphFromSQLite(
    organism: str = "hsa",
    filter_path: Optional[List[str]] = None,
    sqlite_path: Optional[str] = None,
) -> ig.Graph:
    """
    从本地 SQLite 数据库构建经过整理的 KEGG 图。

    Parameters
    ----------
    organism : str
        物种 KEGG 代码，如 ``"hsa"``、``"dre"``、``"mmu"``。
    filter_path : list of str or None
        正则表达式列表，用于过滤掉不需要的通路（如 ``["01100"]`` 过滤总览代谢通路）。
    sqlite_path : str
        SQLite 数据库文件路径。数据库中必须包含 ``item`` 表。

    Returns
    -------
    igraph.Graph
        与 ``buildGraphFromKEGGREST`` 输出结构完全一致的有向图，
        包含顶点属性 ``com``、``NAME``、``entrez`` 以及边属性 ``weight``。
    """
    if sqlite_path is None:
        raise ValueError("必须提供 sqlite_path 参数指向本地 KEGG SQLite 数据库。")
    if not os.path.isfile(sqlite_path):
        raise FileNotFoundError(f"SQLite 数据库未找到: {sqlite_path}")

    categories = ["pathway", "module", "enzyme", "reaction", "compound"]
    conn = sqlite3.connect(sqlite_path)
    cursor = conn.cursor()

    # 校验 item 表存在
    cursor.execute(
        "SELECT name FROM sqlite_master WHERE type='table' AND name='item'"
    )
    if not cursor.fetchone():
        conn.close()
        raise ValueError(f"数据库 {sqlite_path} 中缺少必需的 'item' 表。")

    print(f"Building through local SQLite: {sqlite_path}")

    # ------------------------------------------------------------------
    # 1. 获取 5 类节点的 id → name 映射
    # ------------------------------------------------------------------
    print("Step 1/4: Loading KEGG entries from SQLite...")
    list_list: Dict[str, Dict[str, str]] = {}
    for category in categories:
        print(f"  -> list/{category}...")
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

    # id → category 映射
    map_category = {}
    for category in categories:
        for entry_id in list_list[category]:
            map_category[entry_id] = category

    # ------------------------------------------------------------------
    # 2. 直接链接（边）
    # ------------------------------------------------------------------
    print("Step 2/4: Extracting direct links from SQLite...")
    link_frames = []
    for j in range(len(categories)):
        for i in range(len(categories)):
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
    # 3. 基因层面的映射（用于推断 pathway/module → enzyme）
    # ------------------------------------------------------------------
    print("Step 3/4: Extracting gene-centric mappings from SQLite...")
    m_path_gene, m_mod_gene, m_gene_enzyme, kegg_gene2entrez = _sqlite_gene_mappings(
        cursor, organism
    )

    # enzyme → gene（反向映射）
    m_enzyme_gene: Dict[str, List[str]] = defaultdict(list)
    for gene, enzymes in m_gene_enzyme.items():
        for e in enzymes:
            m_enzyme_gene[e].append(gene)
    m_enzyme_gene = dict(m_enzyme_gene)

    # enzyme → entrez（唯一且排序）
    m_enzyme_gene_entrez = {}
    for enz, genes in m_enzyme_gene.items():
        ncbis = set()
        for g in genes:
            for n in kegg_gene2entrez.get(g, []):
                ncbis.add(n)
        m_enzyme_gene_entrez[enz] = sorted(list(ncbis))

    print("  Done.")

    # ------------------------------------------------------------------
    # 4. 推断连接（pathway / module → enzyme）
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
    # 5. 组装边表（以下逻辑与 REST 版完全一致）
    # ------------------------------------------------------------------
    noinfere_frames = []
    for df in link_frames:
        a_from = df.attrs.get("from")
        a_to = df.attrs.get("to")
        if a_from == "enzyme" and a_to in ("module", "pathway"):
            continue
        noinfere_frames.append(df)

    df_noinfere = (
        pd.concat(noinfere_frames, ignore_index=True)
        if noinfere_frames
        else pd.DataFrame(columns=["from", "to"])
    )
    df_infere = (
        pd.concat(con_infere, ignore_index=True)
        if con_infere
        else pd.DataFrame(columns=["from", "to"])
    )

    df_edges = pd.concat([df_noinfere, df_infere], ignore_index=True)
    df_edges = df_edges.dropna(subset=["from", "to"]).reset_index(drop=True)

    print("Done.")

    # ------------------------------------------------------------------
    # 6. 建图（以下逻辑与 REST 版完全一致）
    # ------------------------------------------------------------------
    print("Building graph...")
    edge_tuples = list(zip(df_edges["from"], df_edges["to"]))
    if not edge_tuples:
        conn.close()
        raise RuntimeError("从 SQLite 中未提取到任何边，无法建图。")

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

    infere_from = (
        set(df_infere["from"].dropna().unique()) if not df_infere.empty else set()
    )
    bad_enzymes = [
        i
        for i, (c, n) in enumerate(zip(g_raw.vs["com"], g_raw.vs["name"]))
        if c == 3 and n not in infere_from
    ]
    if bad_enzymes:
        g_raw.delete_vertices(bad_enzymes)

    org_modules = set(m_mod_gene.keys())
    bad_modules = [
        i
        for i, (c, n) in enumerate(zip(g_raw.vs["com"], g_raw.vs["name"]))
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
        i
        for i, c in enumerate(g_raw.vs["com"])
        if c == 4 and i not in w3_sources
    ]
    if bad_reactions:
        g_raw.delete_vertices(bad_reactions)

    w1_sources = {e.source for e in g_raw.es if e["weight"] == 1}
    bad_compounds = [
        i
        for i, c in enumerate(g_raw.vs["com"])
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
            g_raw.delete_vertices(
                [v.index for v in g_raw.vs if v["name"] in names_out]
            )

    g_raw = _largestcc(g_raw)
    print("Done.")

    # ------------------------------------------------------------------
    # 7. 剪枝冗余传递边（以下逻辑与 REST 版完全一致）
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
    g_curated = g_raw.subgraph_edges(
        edges_by_weight[first_w], delete_vertices=False
    )

    for w in sorted_weights[1:]:
        print(f"Current weight: {w} out of 4...")
        dist_matrix = g_curated.distances(mode="out")
        list_edges = edges_by_weight[w]
        list_ends = [
            (g_raw.es[eid].source, g_raw.es[eid].target) for eid in list_edges
        ]

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
    # 8. 附加顶点元数据
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
        f"SQLite backend; organism={organism}; built from {sqlite_path}"
    )

    conn.close()
    print("Done.")
    return g_curated
