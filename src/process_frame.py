#!/usr/bin/env python
#-------------------------------------------------------------------
# This module is part of coor2graph
# Process individual adjacency matrix
#-------------------------------------------------------------------
#> author: Asdrubal Lozada
#  update: 4/11/2026
#  license: GPL-v3 

import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.ticker as mticker
import sys, os


def draw(G, pos, measures, output_file, shape_intervals=None):
    plt.figure(figsize=(12, 12))

    nodelist = list(measures.keys())
    values = np.array(list(measures.values()))

    sizes = np.full(len(values), 100)

    norm = mcolors.Normalize(
        vmin=values.min(),
        vmax=values.max()
    )

    cmap = plt.cm.RdBu_r


    if shape_intervals is None:
        shape_intervals = []

    assigned_nodes = set()

    colorbar_nodes = None


    for start, end, shape in shape_intervals:

        selected_nodes = [
            node for node in nodelist
            if start <= node <= end
        ]

        if not selected_nodes:
            continue

        selected_values = np.array([
            measures[node]
            for node in selected_nodes
        ])

        selected_sizes = np.full(
            len(selected_nodes),
            100
        )

        nodes = nx.draw_networkx_nodes(
            G,
            pos,
            nodelist=selected_nodes,
            node_size=selected_sizes,
            cmap=cmap,
            node_color=selected_values,
            node_shape=shape,
            edgecolors='black',
            linewidths=0.5,
            alpha=0.9
        )

        nodes.set_norm(norm)

        if colorbar_nodes is None:
            colorbar_nodes = nodes

        assigned_nodes.update(selected_nodes)

    remaining_nodes = [
        node for node in nodelist
        if node not in assigned_nodes
    ]

    if remaining_nodes:

        remaining_values = np.array([
            measures[node]
            for node in remaining_nodes
        ])

        remaining_sizes = np.full(
            len(remaining_nodes),
            100
        )

        nodes = nx.draw_networkx_nodes(
            G,
            pos,
            nodelist=remaining_nodes,
            node_size=remaining_sizes,
            cmap=cmap,
            node_color=remaining_values,
            node_shape='o',
            edgecolors='black',
            linewidths=0.5,
            alpha=0.9
        )

        nodes.set_norm(norm)

        if colorbar_nodes is None:
            colorbar_nodes = nodes

    nx.draw_networkx_edges(
        G,
        pos,
        edge_color="black"
    )


    if colorbar_nodes is not None:

        cbar = plt.colorbar(colorbar_nodes)

        from matplotlib.ticker import MaxNLocator, ScalarFormatter

        cbar.ax.yaxis.set_major_locator(
            MaxNLocator(nbins=10)
        )

        formatter = ScalarFormatter(useMathText=True)
        formatter.set_scientific(True)
        formatter.set_powerlimits((-3, 3))

        cbar.ax.yaxis.set_major_formatter(formatter)

        cbar.ax.tick_params(
            labelsize=30,
            colors="black"
        )


    plt.axis("off")

    plt.subplots_adjust(
        left=0.05,
        right=0.95,
        top=0.95,
        bottom=0.05
    )

    plt.savefig(
        output_file,
        bbox_inches='tight',
        pad_inches=0.1,
        dpi=300
    )

    plt.close()



def process_frame(frame, measure_type='degree', measure_name='Measure', draw_graph=True,shape_intervals=None):
    fstr = str(frame)
    matrix_file = f"adjacency_{fstr}.csv"
    if not os.path.exists(matrix_file):
        print(f"File {matrix_file} not found.")
        return

    print(frame)
    adj_matrix = pd.read_csv(matrix_file, header=None).values
    G = nx.from_numpy_array(adj_matrix)
    pos = nx.spring_layout(G, seed=42)

    if measure_type == 'degree':
        measures = dict(G.degree())
        out_txt = f"degree_{fstr}.txt"
    elif measure_type == 'closeness':
        measures = nx.closeness_centrality(G)
        out_txt = f"closeness_{fstr}.txt"
    elif measure_type == 'betweenness':
        measures = nx.betweenness_centrality(G)
        out_txt = f"betweenness_{fstr}.txt"
    elif measure_type == 'katz':
        measures = nx.katz_centrality(G)
        out_txt = f"katz_{fstr}.txt"
    elif measure_type == 'eigenvector':
        measures = nx.eigenvector_centrality(G)
        out_txt = f"eigenvector_{fstr}.txt"

    with open(out_txt, "w") as f:
        for node, val in measures.items():
            f.write(f"{node}\t{val:.10f}\n")

    if draw_graph:

        out_png = str(measure_type)+'_'+str(fstr)+'.png'
        draw(G, pos, measures, out_png,shape_intervals=shape_intervals)

if __name__ == "__main__":
    process_frame(1)  #

