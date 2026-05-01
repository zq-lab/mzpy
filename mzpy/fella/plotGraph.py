"""
Python port of FELLA::plotGraph

Plots a solution graph from diffusion / pagerank analysis.
When matplotlib is available, it is used as the backend so that
legends are drawn and Jupyter can display the figure inline.
Otherwise falls back to the native igraph cairo plot.
"""

from typing import Optional

import numpy as np
import igraph as ig

from .plotLegend import plotLegend as _plotLegend


def _add_legend_matplotlib(ax, go_annot: bool = False):
    """
    Add R-style legends to a matplotlib Axes.

    .. deprecated::
        This helper is deprecated and will be removed in a future release.
        Use :func:`plotLegend` instead.
    """
    from matplotlib.lines import Line2D

    # Bottom-left: categories
    print(
        "WARNING: _add_legend_matplotlib is deprecated and will be removed. "
        "Use plotLegend instead."
    )

    cat_items = [
        ("Pathway", "#CD0000", "o"),
        ("Module", "#CD96CD", "o"),
        ("Enzyme", "#FFD500", "o"),
        ("Reaction", "#8DB6CD", "o"),
        ("Compound", "#548B54", "o"),
        ("Input compound", "#548B54", "s"),
    ]
    cat_handles = [
        Line2D(
            [0],
            [0],
            marker=shape,
            color="w",
            markerfacecolor=color,
            markeredgecolor="black",
            markersize=12,
            label=label,
        )
        for label, color, shape in cat_items
    ]
    leg1 = ax.legend(
        handles=cat_handles,
        loc="lower left",
        title="Categories for each node",
        ncol=3,
        frameon=True,
        fancybox=True,
        fontsize=11,
        title_fontsize=11,
    )
    ax.add_artist(leg1)

    if go_annot:
        go_items = [
            ("Simil < 0.5", "#FFD500"),
            ("Simil < 0.7", "#FF5500"),
            ("Simil < 0.9", "#FF0000"),
            ("Simil <= 1", "#B300FF"),
        ]
        go_handles = [
            Line2D(
                [0],
                [0],
                marker="^",
                color="w",
                markerfacecolor=color,
                markeredgecolor="black",
                markersize=12,
                label=label,
            )
            for label, color in go_items
        ]
        ax.legend(
            handles=go_handles,
            loc="lower right",
            title="Enzymes with CC similarity",
            ncol=2,
            frameon=True,
            fancybox=True,
            fontsize=10,
            title_fontsize=10,
        )


def plotGraph(
    graph: Optional[ig.Graph] = None,
    layout: bool = False,
    graph_layout: Optional[np.ndarray] = None,
    NamesAsLabels: bool = False,
    target: Optional[str] = None,
    plotLegend: bool = True,
    bbox=(1200, 1000),
    **kwargs,
):
    """
    Plot a solution graph.

    Parameters
    ----------
    graph : igraph.Graph
        Solution graph.
    layout : bool
        If True, return layout DataFrame.
    graph_layout : np.ndarray or None
        Pre-computed 2-D layout matrix.
    NamesAsLabels : bool
        Use KEGG names instead of IDs as labels.
    target : str or None
        Output file path (e.g. "plot.png"). If None, the figure is
        returned for display (e.g. in Jupyter).
    plotLegend : bool
        Whether to draw the category legend (R-style).

    Returns
    -------
    pd.DataFrame or matplotlib.figure.Figure or None
    """
    if graph is None or graph.vcount() == 0:
        print("The graph is empty and won't be plotted.")
        return None

    go_annot = "GO.simil" in graph.vertex_attributes()
    go_simil = [v.attributes().get("GO.simil", None) for v in graph.vs] if go_annot else [None] * graph.vcount()

    coms = [v["com"] for v in graph.vs]
    inputs = [v.attributes().get("input", False) for v in graph.vs]

    if graph_layout is None:
        graph_layout = graph.layout("auto")
    else:
        graph_layout = ig.Layout(graph_layout.tolist() if hasattr(graph_layout, "tolist") else graph_layout)

    coords = np.array(graph_layout.coords)
    if coords[:, 0].max() != coords[:, 0].min():
        coords[:, 0] = (coords[:, 0] - coords[:, 0].min()) / (coords[:, 0].max() - coords[:, 0].min()) * 2 - 1
    if coords[:, 1].max() != coords[:, 1].min():
        coords[:, 1] = (coords[:, 1] - coords[:, 1].min()) / (coords[:, 1].max() - coords[:, 1].min()) * 2 - 1
    graph_layout = ig.Layout(coords.tolist())

    # Colors / shapes / frame colors
    map_solid = {1: "#CD0000", 2: "#CD96CD", 3: "#FFA200", 4: "#8DB6CD", 5: "#548B54"}
    colors = []
    frame_colors = []
    shapes = []
    for i, v in enumerate(graph.vs):
        solid = map_solid.get(coms[i], "grey")
        if go_annot and go_simil[i] is not None:
            gs = go_simil[i]
            if gs < 0.5:
                solid = "#FFD500"
            elif gs < 0.7:
                solid = "#FF5500"
            elif gs < 0.9:
                solid = "#FF0000"
            else:
                solid = "#B300FF"
            frame_colors.append("#CD0000")
            shapes.append("triangle")
        else:
            frame_colors.append("black")
            shapes.append("square" if inputs[i] else "circle")
        colors.append(solid)

    # Sizes
    map_size = {1: 20, 2: 16, 3: 13, 4: 11, 5: 9}
    sizes = [map_size.get(c, 9) for c in coms]
    sizes = [s * (300 / graph.vcount()) ** (1 / 3) * 0.75 for s in sizes]

    labels = [
        v["label"] if NamesAsLabels and "label" in v.attributes() else v["name"]
        for v in graph.vs
    ]

    # Synchronise the vertex label attribute so both matplotlib and cairo
    # backends display the intended text (IDs when NamesAsLabels=False).
    graph.vs["label"] = labels

    visual_style = {
        "layout": graph_layout,
        "vertex_color": colors,
        "vertex_shape": shapes,
        "vertex_size": sizes,
        "vertex_label": labels,
        "vertex_label_size": 12,
        "vertex_label_dist": 0.12 * (300 / graph.vcount()) ** (1 / 3),
        "vertex_label_color": "black",
        "vertex_frame_color": frame_colors,
        "edge_color": "#000000AA",
        "edge_width": 0.5,
        "edge_arrow_size": 0.25,
        "bbox": bbox,
        "margin": 50,
    }
    visual_style.update(kwargs)

    # Try matplotlib backend first (supports legend & Jupyter inline display)
    try:
        import matplotlib.pyplot as plt

        style = visual_style.copy()
        bbox = style.pop("bbox", bbox)
        style.pop("margin", None)
        style.pop("vertex_label_dist", None)

        # Map shapes for matplotlib compatibility
        if "vertex_shape" in style:
            style["vertex_shape"] = [
                "triangle-up" if s == "triangle" else s
                for s in style["vertex_shape"]
            ]

        fig, ax = plt.subplots(figsize=(bbox[0] / 100, bbox[1] / 100))
        ig.plot(graph, target=ax, **style)

        ax.axis("off")

        if plotLegend:
            _plotLegend(GO_annot=go_annot, ax=ax)

        if target is not None:
            fig.savefig(target, bbox_inches="tight", pad_inches=0.1)
            plt.close(fig)
            out = None
        else:
            out = fig

        if layout:
            import pandas as pd

            prefixes = {1: "", 2: "md:", 3: "ec:", 4: "rn:", 5: "cpd:"}
            df = pd.DataFrame({
                "x": coords[:, 0],
                "y": coords[:, 1],
                "out.id": [prefixes.get(c, "") + v["name"] for c, v in zip(coms, graph.vs)],
                "out.name": labels,
            })
            return df
        return out

    except Exception as e:
        # Fall back to native cairo plot
        try:
            import matplotlib.pyplot as plt

            plt.close("all")
        except Exception:
            pass
        print(f"Matplotlib backend failed ({e}), falling back to cairo.")

    # --- cairo fallback ---
    plot = ig.plot(graph, target=target, **visual_style)
    if target is None and hasattr(plot, "show"):
        plot.show()

    if layout:
        import pandas as pd

        prefixes = {1: "", 2: "md:", 3: "ec:", 4: "rn:", 5: "cpd:"}
        out = pd.DataFrame({
            "x": coords[:, 0],
            "y": coords[:, 1],
            "out.id": [prefixes.get(c, "") + v["name"] for c, v in zip(coms, graph.vs)],
            "out.name": labels,
        })
        return out
    return None
