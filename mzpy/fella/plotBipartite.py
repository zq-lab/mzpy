"""
Python port of FELLA:::plotBipartite

Plots a bipartite graph tailored for the hypergeometric test.
When matplotlib is available, it is used for Jupyter inline display.
"""
from typing import Optional

import numpy as np
import igraph as ig


def plotBipartite(
    graph: Optional[ig.Graph] = None,
    layout: bool = False,
    target: Optional[str] = None,
    **kwargs,
):
    """
    Plot a bipartite graph.

    Parameters
    ----------
    graph : igraph.Graph
        Bipartite graph from hypergeometric test.
    layout : bool
        If True, return layout DataFrame.
    target : str or None
        Output file path (e.g. "plot.png"). If None, the figure is
        returned for display (e.g. in Jupyter).

    Returns
    -------
    pd.DataFrame or matplotlib.figure.Figure or None
    """
    if graph is None or graph.vcount() == 0:
        print("Graph is empty.")
        return None

    coms = [v["com"] for v in graph.vs]
    types = [c == 5 for c in coms]

    # Bipartite layout
    graph_layout = graph.layout_bipartite(types=types, hgap=0.5, vgap=2.5)
    coords = np.array(graph_layout.coords)
    # Swap x/y to match original orientation
    coords = coords[:, [1, 0]]
    # Normalise roughly to -1..1
    if coords[:, 0].max() != coords[:, 0].min():
        coords[:, 0] = (coords[:, 0] - coords[:, 0].min()) / (coords[:, 0].max() - coords[:, 0].min()) * 2 - 1
    if coords[:, 1].max() != coords[:, 1].min():
        coords[:, 1] = (coords[:, 1] - coords[:, 1].min()) / (coords[:, 1].max() - coords[:, 1].min()) * 2 - 1
    graph_layout = ig.Layout(coords.tolist())

    colors = ["palegreen4" if c == 5 else "red3" for c in coms]
    shapes = ["square" if c == 5 else "circle" for c in coms]
    labels = [v["name"] for v in graph.vs]
    label_colors = colors

    visual_style = {
        "layout": graph_layout,
        "vertex_color": colors,
        "vertex_shape": shapes,
        "vertex_size": 12,
        "vertex_label": labels,
        "vertex_label_dist": 0.8,
        "vertex_label_color": label_colors,
        "edge_color": "#000000AA",
        "edge_width": 0.5,
        "edge_arrow_size": 0.25,
        "bbox": (700, 500),
        "margin": 60,
    }
    visual_style.update(kwargs)

    # Try matplotlib backend for Jupyter display
    try:
        import matplotlib.pyplot as plt

        style = visual_style.copy()
        bbox = style.pop("bbox", (700, 500))
        style.pop("margin", None)
        style.pop("vertex_label_dist", None)

        fig, ax = plt.subplots(figsize=(bbox[0] / 100, bbox[1] / 100))
        ig.plot(graph, target=ax, **style)
        ax.axis("off")

        if target is not None:
            fig.savefig(target, bbox_inches="tight", pad_inches=0.1)
            plt.close(fig)
            out = None
        else:
            out = fig

        if layout:
            import pandas as pd
            out_df = pd.DataFrame({
                "x": coords[:, 0],
                "y": coords[:, 1],
                "out.id": [
                    f"cpd:{v['name']}" if v["com"] == 5 else v["name"]
                    for v in graph.vs
                ],
            })
            return out_df
        return out

    except Exception as e:
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
        out = pd.DataFrame({
            "x": coords[:, 0],
            "y": coords[:, 1],
            "out.id": [
                f"cpd:{v['name']}" if v["com"] == 5 else v["name"]
                for v in graph.vs
            ],
        })
        return out
    return None
