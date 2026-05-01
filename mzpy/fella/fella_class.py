"""
High-level wrapper class for FELLA enrichment analysis.

This class provides a user-friendly interface to the fella-py package,
encapsulating the typical workflow of loading a database, running enrichment,
and extracting / plotting / exporting results.

All threshold / filtering parameters are fixed at ``enrich()`` time so that
subsequent result-extraction methods (table, graph, enzymes, plot) always
operate on the same settings and produce consistent outputs.
"""

from typing import List, Optional, Any, Union, Dict
import os

import pandas as pd
import igraph as ig

from .loadKEGGdata import loadKEGGdata
from .enrich import enrich as _enrich
from .generateResultsGraph import generateResultsGraph
from .generateResultsTable import generateResultsTable
from .generateEnzymesTable import generateEnzymesTable
from .exportResults import exportResults
from .AllMethods import show, plot
from .get_ import getGraph, getInput, getPscores, getValid, getStatus
from .list_ import listMethods


class Fella:
    """
    High-level interface for FELLA metabolomics enrichment analysis.

    Parameters
    ----------
    organism : str
        KEGG organism code (e.g. ``'hsa'``, ``'mmu'``, ``'dre'``).
    h5_path : str
        Path to the HDF5 database file produced by ``buildDataFromGraph``.
    load_matrix : list of str or None, optional
        Which heavy matrices to load into memory.
        Can contain ``"diffusion"`` and/or ``"pagerank"``.
        If None, no heavy matrices are loaded.

    Examples
    --------
    >>> fella = Fella(organism="dre", h5_path="./data/dre.h5",
    ...               load_matrix=["diffusion", "pagerank"])
    >>> results = fella.enrich(
    ...     compounds=["C00002", "C00040", "C00186"],
    ...     methods="diffusion",
    ...     threshold=0.05,
    ...     output_dir="./results"
    ... )
    """

    # ------------------------------------------------------------------
    # Initialisation
    # ------------------------------------------------------------------

    def __init__(
        self,
        organism: str,
        h5_path: str,
        load_matrix: Optional[List[str]] = None,
    ):
        self.organism = organism
        self.h5_path = h5_path
        self.data = loadKEGGdata(h5_path, load_matrix=load_matrix)

        # Enrichment result (must be computed first)
        self.obj = None
        self._last_method: Optional[str] = None
        self._methods_ran: set = set()

        # Threshold / filter parameters fixed at enrich() time
        self._threshold = 0.05
        self._nlimit = 250
        self._plimit = 15
        self._thresholdConnectedComponent = 0.05
        self._minDegree = 0

        # Cached intermediate results – all start as None
        self._results_graph = None
        self._results_graph_key = None

        self._results_table = None
        self._results_table_key = None

        self._enzymes_table = None
        self._enzymes_table_key = None

    # ------------------------------------------------------------------
    # Cache helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _make_key(method: str, kwargs: dict):
        """Build an immutable cache key from method name and kwargs."""
        return (method, tuple(sorted(kwargs.items())))

    def _clear_caches(self):
        """Invalidate all cached intermediate results."""
        self._results_graph = None
        self._results_graph_key = None
        self._results_table = None
        self._results_table_key = None
        self._enzymes_table = None
        self._enzymes_table_key = None

    # ------------------------------------------------------------------
    # Core workflow
    # ------------------------------------------------------------------

    def enrich(
        self,
        compounds: List[str],
        compoundsBackground: Optional[List[str]] = None,
        methods: Optional[Union[str, List[str]]] = None,
        threshold: float = 0.05,
        nlimit: int = 250,
        plimit: int = 15,
        thresholdConnectedComponent: float = 0.05,
        minDegree: int = 0,
        output_dir: Optional[str] = None,
        **kwargs,
    ) -> Dict[str, Dict[str, Any]]:
        """
        Run enrichment analysis and generate all result artefacts.

        All filtering parameters (threshold, nlimit, etc.) are fixed at this
        call.  Subsequent ``generate_results_*`` methods and ``plot()`` reuse
        these values automatically so that table, graph and enzymes are always
        consistent.

        If *output_dir* is provided, three files are written per method:

        - ``{method}_results_table.csv``
        - ``{method}_enzymes_table.csv`` (diffusion / pagerank only)
        - ``{method}_network.svg``

        Parameters
        ----------
        compounds : list of str
            KEGG compound IDs to enrich.
        compoundsBackground : list of str or None
            Custom background compound IDs.
        methods : str or list of str or None
            Enrichment method(s). Supported: ``"hypergeom"``, ``"diffusion"``,
            ``"pagerank"``. If None, all available methods are run.
        threshold : float
            p-score threshold for result inclusion (default 0.05).
        nlimit : int
            Max nodes for diffusion / pagerank (default 250).
        plimit : int
            Max pathways for hypergeom (default 15).
        thresholdConnectedComponent : float
            Minimum relative size of connected components to keep in the
            network graph (default 0.05).
        minDegree : int
            Minimum vertex degree for the result graph (default 0).
        output_dir : str or None
            If given, CSV tables and SVG network plots are saved here.
        **kwargs
            Additional arguments passed to the underlying enrichment runners
            (e.g. ``approx``, ``niter``, ``t_df``).

        Returns
        -------
        dict
            Nested dictionary ``{method: {"table": df, "enzymes": df|None,
            "graph": ig.Graph, "figure": matplotlib.Figure|None}}``.
        """
        if isinstance(methods, str):
            methods = [methods]

        # New enrichment run -> previous caches are stale
        print('clear caches ... ')
        self._clear_caches()

        # Fix filtering parameters for this enrichment session
        self._threshold = threshold
        self._nlimit = nlimit
        self._plimit = plimit
        self._thresholdConnectedComponent = thresholdConnectedComponent
        self._minDegree = minDegree
        
        self.obj = _enrich(
            compounds=compounds,
            compoundsBackground=compoundsBackground,
            methods=methods,
            data=self.data,
            **kwargs,
        )

        if methods is not None:
            self._methods_ran.update(methods)
            self._last_method = methods[0]
        else:
            self._methods_ran.update(listMethods())
            self._last_method = "diffusion"

        # Generate results for every method that was run
        results: Dict[str, Dict[str, Any]] = {}
        for m in sorted(self._methods_ran):
            res: Dict[str, Any] = {}
            res["table"] = self.generate_results_table(method=m)
            if m in ("diffusion", "pagerank"):
                res["enzymes"] = self.generate_enzymes_table(method=m)
            else:
                res["enzymes"] = None
            res["graph"] = self.generate_results_graph(method=m)
            res["figure"] = None
            if res["graph"] is not None and res["graph"].vcount() > 0:
                res["figure"] = self.plot_graph(graph=res["graph"], method=m)
            results[m] = res

        # Save to disk if requested
        if output_dir is not None:
            os.makedirs(output_dir, exist_ok=True)
            for m, res in results.items():
                prefix = os.path.join(output_dir, m)
                if res["table"] is not None:
                    res["table"].to_csv(f"{prefix}_results_table.csv", index=False)
                if res["enzymes"] is not None:
                    res["enzymes"].to_csv(f"{prefix}_enzymes_table.csv", index=False)
                if res["figure"] is not None:
                    res["figure"].savefig(
                        f"{prefix}_network.svg",
                        bbox_inches="tight",
                        pad_inches=0.1,
                    )
                    import matplotlib.pyplot as plt
                    plt.close(res["figure"])
            print(f"Results saved to {output_dir}")

        return results

    # ------------------------------------------------------------------
    # Result extraction  (cached, parameter-free)
    # ------------------------------------------------------------------

    def generate_results_graph(
        self,
        method: Optional[str] = None,
    ) -> Optional[ig.Graph]:
        """
        Generate a sub-network graph with the most significant nodes.

        Uses the filtering parameters fixed by the last ``enrich()`` call.

        Parameters
        ----------
        method : str or None
            Enrichment method. If None, the last run method is used.

        Returns
        -------
        igraph.Graph or None
        """
        self._check_enriched()
        m = self._resolve_method(method)
        params = {
            "threshold": self._threshold,
            "nlimit": self._nlimit,
            "plimit": self._plimit,
            "thresholdConnectedComponent": self._thresholdConnectedComponent,
            "minDegree": self._minDegree,
        }
        key = self._make_key(m, params)

        if self._results_graph is not None and self._results_graph_key == key:
            return self._results_graph

        graph = generateResultsGraph(
            object=self.obj, method=m, data=self.data, **params
        )

        if graph is not None and graph.vcount() > 0:
            self._results_graph = graph
            self._results_graph_key = key
        return graph

    def generate_results_table(
        self,
        method: Optional[str] = None,
    ) -> Optional[pd.DataFrame]:
        """
        Generate a results table from the enrichment analysis.

        Uses the *threshold* fixed by the last ``enrich()`` call.

        Parameters
        ----------
        method : str or None
            Enrichment method. If None, the last run method is used.

        Returns
        -------
        pd.DataFrame or None
        """
        self._check_enriched()
        m = self._resolve_method(method)
        params = {
            "threshold": self._threshold,
            "nlimit": self._nlimit,
            "plimit": self._plimit,
        }
        key = self._make_key(m, params)

        if self._results_table is not None and self._results_table_key == key:
            return self._results_table

        df = generateResultsTable(
            object=self.obj, method=m, data=self.data, **params
        )

        # Custom sort: Entry type by biological hierarchy, then p.score asc.
        if df is not None and "Entry type" in df.columns:
            desired_order = ["pathway", "module", "reaction", "enzyme", "compound"]
            present = [c for c in desired_order if c in df["Entry type"].values]
            df["Entry type"] = pd.Categorical(
                df["Entry type"], categories=present, ordered=True
            )
            df = df.sort_values(by=["Entry type", "p.score"]).reset_index(drop=True)

        self._results_table = df
        self._results_table_key = key
        return self._results_table

    def generate_enzymes_table(
        self,
        method: Optional[str] = None,
    ) -> Optional[pd.DataFrame]:
        """
        Generate a table of enzyme families and their annotated genes.

        Uses the *threshold* fixed by the last ``enrich()`` call.

        Parameters
        ----------
        method : str or None
            Enrichment method. Must be ``"diffusion"`` or ``"pagerank"``.
            If None, the last run method is used.

        Returns
        -------
        pd.DataFrame or None
        """
        self._check_enriched()
        m = self._resolve_method(method)
        if m not in ("diffusion", "pagerank"):
            raise ValueError(
                f"Enzymes table is only available for 'diffusion' or 'pagerank', "
                f"but got '{m}'."
            )
        params = {
            "threshold": self._threshold,
            "nlimit": self._nlimit,
        }
        key = self._make_key(m, params)

        if self._enzymes_table is not None and self._enzymes_table_key == key:
            return self._enzymes_table

        self._enzymes_table = generateEnzymesTable(
            object=self.obj, method=m, data=self.data, **params
        )
        self._enzymes_table_key = key
        return self._enzymes_table

    # ------------------------------------------------------------------
    # Plotting
    # ------------------------------------------------------------------

    def plot(
        self,
        method: Optional[str] = None,
        layout: bool = False,
        **kwargs,
    ) -> Any:
        """
        Plot the enrichment results.

        Parameters
        ----------
        method : str or None
            Enrichment method. If None, the last run method is used.
        layout : bool
            If True, return layout DataFrame.
        **kwargs
            Plotting-only arguments passed to the backend
            (e.g. ``target``, ``plotLegend``, ``NamesAsLabels``).
            ``useNamesAsLabels`` is accepted as an alias for
            ``NamesAsLabels`` for backward compatibility.

        Returns
        -------
        matplotlib.figure.Figure, pd.DataFrame or None
        """
        self._check_enriched()
        m = self._resolve_method(method)

        plot_kwargs = dict(kwargs)
        plot_kwargs["layout"] = layout

        # Alias for backward compatibility
        if "useNamesAsLabels" in plot_kwargs:
            plot_kwargs["NamesAsLabels"] = plot_kwargs.pop("useNamesAsLabels")

        graph = self.generate_results_graph(method=m)
        if graph is None or graph.vcount() == 0:
            print("Graph is empty, nothing to plot.")
            return None

        if m == "hypergeom":
            from .plotBipartite import plotBipartite
            return plotBipartite(graph=graph, **plot_kwargs)
        else:
            from .plotGraph import plotGraph
            return plotGraph(graph=graph, **plot_kwargs)

    def plot_graph(
        self,
        graph: Optional[ig.Graph] = None,
        method: Optional[str] = None,
        **kwargs,
    ) -> Any:
        """
        Plot a solution graph.

        If *graph* is not provided, a result graph is first generated from
        the current enrichment using *method*.

        Parameters
        ----------
        graph : igraph.Graph or None
            Graph to plot. If None, ``generate_results_graph()`` is called.
        method : str or None
            Enrichment method used when *graph* is None.
        **kwargs
            Additional arguments passed to ``plotGraph``.

        Returns
        -------
        matplotlib.figure.Figure, pd.DataFrame or None
        """
        from .plotGraph import plotGraph

        if graph is None:
            graph = self.generate_results_graph(method=method)

        if graph is None or graph.vcount() == 0:
            print("Graph is empty, nothing to plot.")
            return None
        return plotGraph(graph=graph, **kwargs)

    # ------------------------------------------------------------------
    # Export
    # ------------------------------------------------------------------

    def export_results(
        self,
        file: str,
        format: str = "csv",
        method: Optional[str] = None,
    ):
        """
        Export enrichment results to a file.

        Parameters
        ----------
        file : str
            Output file path.
        format : str
            Export format: ``"csv"``, ``"enzyme"``, ``"igraph"``, or any
            igraph-supported graph format.
        method : str or None
            Enrichment method. If None, the last run method is used.
        """
        self._check_enriched()
        m = self._resolve_method(method)

        if format == "csv":
            df = self.generate_results_table(method=m)
            if df is not None:
                df.to_csv(file, index=False)
        elif format == "enzyme":
            if m not in ("diffusion", "pagerank"):
                raise ValueError(
                    f"Enzymes are only reported in diffusion and pagerank, "
                    f"but not in method {m}"
                )
            df = self.generate_enzymes_table(method=m)
            if df is not None:
                df.to_csv(file, index=False)
        else:
            graph = self.generate_results_graph(method=m)
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

    # ------------------------------------------------------------------
    # Utilities / introspection
    # ------------------------------------------------------------------

    def show(self):
        """Print a summary of the loaded database or enrichment results."""
        if self.obj is None:
            show(self.data)
        else:
            show(self.obj)

    def summary(self) -> dict:
        """Return a dict summarising the current state."""
        info = {
            "organism": self.organism,
            "h5_path": self.h5_path,
            "data_loaded": getStatus(self.data) == "loaded",
            "graph_nodes": 0,
            "graph_edges": 0,
            "enrichment_run": self.obj is not None,
            "methods_ran": sorted(self._methods_ran),
            "last_method": self._last_method,
            "cached_results_graph": self._results_graph is not None,
            "cached_results_table": self._results_table is not None,
            "cached_enzymes_table": self._enzymes_table is not None,
            "threshold": self._threshold,
            "nlimit": self._nlimit,
            "plimit": self._plimit,
            "thresholdConnectedComponent": self._thresholdConnectedComponent,
            "minDegree": self._minDegree,
        }
        g = getGraph(self.data)
        if g is not None:
            info["graph_nodes"] = g.vcount()
            info["graph_edges"] = g.ecount()
        if self.obj is not None:
            info["input_compounds"] = getInput(self.obj)
        return info

    def get_graph(self) -> Optional[ig.Graph]:
        """Return the underlying KEGG graph."""
        return getGraph(self.data)

    def get_input(self) -> List[str]:
        """Return the list of input compounds used in the last enrichment."""
        self._check_enriched()
        return getInput(self.obj)

    def get_pscores(self, method: Optional[str] = None):
        """
        Return p-scores / p-values from the enrichment analysis.

        Parameters
        ----------
        method : str or None
            Enrichment method. If None, the last run method is used.

        Returns
        -------
        dict, pd.Series or None
        """
        self._check_enriched()
        m = self._resolve_method(method)
        if not getValid(self.obj, m):
            raise ValueError(f"Method '{m}' has not been executed or failed.")
        return getPscores(self.obj, m)

    def __repr__(self) -> str:
        return (
            f"Fella(organism='{self.organism}', h5_path='{self.h5_path}', "
            f"enriched={self.obj is not None})"
        )

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _check_enriched(self):
        if self.obj is None:
            raise RuntimeError(
                "No enrichment results available. Please call enrich() first."
            )

    def _resolve_method(self, method: Optional[str]) -> str:
        if method is not None:
            return method
        if self._last_method is not None:
            return self._last_method
        raise ValueError(
            "No method specified and no enrichment has been run yet."
        )


# ------------------------------------------------------------------------------
# Example usage
# ------------------------------------------------------------------------------
if __name__ == "__main__":
    """
    示例 1：一站式富集分析 —— enrich() 中固定所有阈值参数并自动保存结果
    ------------------------------------------------------------------------
    所有过滤参数（threshold、minDegree 等）在 enrich() 中一次性设定，
    后续 generate_results_table / generate_results_graph / generate_enzymes_table
    会自动复用这些值，确保表格、网络图和酶基因列表完全一致。
    """

    # 1. 初始化（加载 KEGG 数据库）
    # fella = Fella(
    #     organism="cel",
    #     h5_path="../../../../db/FELLA_py_db/cel.h5",
    #     load_matrix=["diffusion", "pagerank"]
    # )

    # 2. 运行富集分析，固定阈值，同时输出到文件夹
    # results = fella.enrich(
    #     compounds=["C00002", "C00040", "C00186"],
    #     methods="diffusion",
    #     threshold=0.05,
    #     thresholdConnectedComponent=0.3,
    #     minDegree=2,
    #     output_dir="./fella_results"
    # )
    # 自动保存：
    #   ./fella_results/diffusion_results_table.csv
    #   ./fella_results/diffusion_enzymes_table.csv
    #   ./fella_results/diffusion_network.svg

    # 3. 手动访问返回结果（不依赖 output_dir）
    # diffusion_res = results["diffusion"]
    # print(diffusion_res["table"].head())      # 富集结果表
    # print(diffusion_res["enzymes"].head())    # 酶基因表
    # print(diffusion_res["graph"].vcount())    # 网络图节点数
    # diffusion_res["figure"].savefig("network.png")  # 手动保存绘图

    # --------------------------------------------------------------------------
    # 示例 2：分步调用（参数已在 enrich() 中固定，后续方法无需再传）
    # --------------------------------------------------------------------------
    # 先运行富集（不设 output_dir）
    # results = fella.enrich(
    #     compounds=["C00002", "C00040", "C00186"],
    #     methods="diffusion",
    #     threshold=0.05,
    #     minDegree=2
    # )
    #
    # 之后随时复用同一套参数获取结果
    # table = fella.generate_results_table()          # 自动使用 threshold=0.05
    # graph = fella.generate_results_graph()          # 自动使用同样的阈值
    # enz   = fella.generate_enzymes_table()          # 自动使用同样的阈值
    # fig   = fella.plot(target="network.svg")        # 绘图也基于同一套参数
