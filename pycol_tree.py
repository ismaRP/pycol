import json
import networkx as nx
from networkx.readwrite import json_graph
import requests
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm


def get_ancestors(tree, taxid, loa):
    pred = list(tree.predecessors(taxid))
    if not pred:
        return loa
    else:
        loa.append(pred[0])
        get_ancestors(tree, pred[0], loa)
        return loa


def check_ancestor(taxid, anc_check, tree):
    if taxid == anc_check:
        return True
    pred = list(tree.predecessors(taxid))
    if not pred:
        return False
    if anc_check == pred[0]:
        return True
    else:
        return check_ancestor(pred[0], anc_check, tree)


def get_lca(tree, taxid_list):
    if len(taxid_list) == 1:
        return taxid_list[0]
    lca = nx.lowest_common_ancestor(tree, taxid_list[0], taxid_list[1])
    if len(taxid_list) > 2:
        for taxid in taxid_list[2:]:
            lca = nx.lowest_common_ancestor(tree, taxid, lca)
    return lca


def get_tree(taxid_list, path=None):
    if path is None:
        url_base = 'http://api.unipept.ugent.be/api/v1/taxa2tree.json?'
        q = ''.join(['input[]=' + str(t) + '&' for t in taxid_list])
        r = requests.get(url_base + q)
        data = r.json()
        tree = json_graph.tree_graph(data)
    else:
        with open(path, 'r') as f:
            data = json.load(f)
        tree = json_graph.tree_graph(data)

    # Get global lca
    lca = get_lca(tree, taxid_list)

    # Prune tree above global lca
    ancestors = get_ancestors(tree, lca, [])
    tree.remove_nodes_from(ancestors)
    return tree


def export_legend(legend, filename="legend.png"):
    fig = legend.figure
    fig.canvas.draw()
    fig.set_size_inches(10.5, 18.5)
    bbox = legend.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    fig.savefig(filename, dpi=300, bbox_inches=bbox)


def draw_tree_fig(tree, data_label, data_color, log_node, path):
    node_color = np.zeros(len(tree))
    for i, (_, n) in enumerate(tree.nodes(data=data_color)):
        if log_node:
            node_color[i] = np.log10(n) if n > 0 else 0
        else:
            node_color[i] = n

    node_cmap = cm.get_cmap('YlOrRd')
    nodenorm = plt.Normalize(np.min(node_color), np.max(node_color))
    nodesm = cm.ScalarMappable(norm=nodenorm, cmap=node_cmap)
    nodesm.set_array([])

    # Get node labels
    node_labels = {}
    for taxid, name in tree.nodes.data(data_label):
        node_labels[taxid] = name

    fig = plt.figure(figsize=(8, 5))
    ax = fig.add_subplot()
    pos = nx.nx_agraph.pygraphviz_layout(tree, prog='dot')
    label_pos = {}
    for n, p in pos.items():
        label_pos[n] = np.array(p) - 20

    nx.draw_networkx_edges(tree, pos=pos, arrows=False, alpha=0.5, node_size=100)
    nx.draw_networkx_nodes(
        tree, pos=pos, node_color=node_color, node_size=100, alpha=0.8, node_shape='o', ax=ax,
        vmin=np.min(node_color), vmax=np.max(node_color), cmap=node_cmap
    )
    nx.draw_networkx_labels(tree, label_pos, node_labels, ax=ax)

    fig.set_tight_layout(True)
    ax.margins(0.2)
    ax.axis("off")

    cmfig, cmax = plt.subplots()
    cbar = plt.colorbar(nodesm, ax=cmax, orientation='vertical', fraction=0.05, pad=0.05)
    cmax.remove()
    if log_node:
        cbar.ax.set_title('log10(# peptides)')
    else:
        cbar.ax.set_title('# peptides')

    if path is not None:
        fig.savefig(path, dpi=300)
        cmpath = path.split('/')
        cmpath = cmpath[0:len(cmpath)-1]
        cmpath = '/'.join(cmpath) + '/cmap.png'
        cmfig.savefig(cmpath, bbox_inches='tight')
        plt.close()
    else:
        plt.close()
        return fig, cmfig


def save_tree(tree, path):
    # Get root as node with in_degree = 0
    for nid, d in tree.in_degree:
        if d == 0:
            root = nid
            break
    tree_json = nx.tree_data(tree, root)
    with open(path, 'w') as f:
        json.dump(tree_json, f)

