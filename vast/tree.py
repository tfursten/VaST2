from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from plotly.colors import hex_to_rgb, label_rgb
import numpy as np



try:
    from vast.patterns import (
        get_patterns)
    from vast.utils import (
        get_metadata_from_genome
    )
except ImportError:
    from patterns import (
        get_patterns)
    from utils import (
        get_metadata_from_genome
    )


def get_snp_alignment(matrix):
    """
    Given a SNP dataframe with genome columns, create
    an alignment object.
    """
    seqs = []
    for c in matrix:
        seqs.append(SeqRecord(Seq("".join(matrix[c].values)), id=str(c)))
    return MultipleSeqAlignment(seqs)
    

def get_distance_matrix(alignment):
    calculator = DistanceCalculator('identity')
    distance_matrix = calculator.get_distance(alignment)
    return distance_matrix

def get_tree(distance_matrix):
    constructor = DistanceTreeConstructor()
    tree = constructor.nj(distance_matrix)
    tree.rooted = True
    return tree


def vast_target_tree(target_results, metadata=None):
    cols = [c for c in target_results.columns.values]
    if metadata is not None:
        # get metadata labels
        metalabels = [get_metadata_from_genome(metadata, c) for c in cols]
        # only include them if they are not the same as the column labels (i.e for full resolution)
        cols = ["{0} {1}".format(c, m) if c != m else c for c, m in zip(cols, metalabels) ]
    snps = pd.DataFrame(
        get_patterns(target_results.sort_index(), 1, 1)['patterns'],
            columns=cols, dtype=str)
    return get_tree(get_distance_matrix(get_snp_alignment(snps)))


# def get_sankey_data(patterns, metadata=None, alpha=0.5):
#     labels = []
#     sources = []
#     targets = []
#     values = []
#     colors = []

#     if metadata is not None:
#         n_cols = len(metadata['names'])
#         cols = px.colors.sample_colorscale("Portland", [n/(n_cols -1) for n in range(n_cols)])
#         palette = {n: c for n, c in zip(metadata['names'].index, cols)}
#         # Map colors to individual genomes based on metadata group membership
#         cm = {i: palette[g] for i, g in enumerate(metadata['values'].values)}

#     else:
#         cm = {}

#     for i, t in enumerate(patterns):
#         labels += ["{0}-{1}".format(i, v) for v in np.unique(t)] 
            
#     for g, genome in enumerate(patterns.T):
#         for target in range(1, len(genome)):
#             sources.append(
#                 labels.index("{0}-{1}".format(target - 1, genome[target-1])))
#             targets.append(
#                 labels.index("{0}-{1}".format(target, genome[target])))
#             colors.append(cm.get(g, '#bababa'))
#             values.append(1)
    
#     if metadata is not None:
#         # add genome labels to the last target
#         last_group = np.unique(patterns[-1], return_inverse=True)[1]
        
#         for i in range(len(labels)):
#             l = labels[i]
#             print(l, len(patterns))
#             if int(l.split("-")[0]) < len(patterns) - 1:
#                 labels[i] = ""
#             else:
#                 group = int(l.split("-")[1])
#                 print(group)
#                 locs = np.where(last_group == group)[0]
#                 names = metadata['values'].index.values[locs]
#                 labels[i] = "<br>".join(names)
    

#     print(np.unique(patterns[-1], return_inverse=True))
#     return {
#         'labels': labels,
#         'sources': sources,
#         'targets': targets,
#         'values': values,
#         'colors': colors
#     }      


# def draw_sankey_diagram(patterns, outfile, metadata=None):
#     sankey_data = get_sankey_data(patterns, metadata)
#     print(sankey_data['labels'])
#     fig = go.Figure(data=[go.Sankey(
#         node = dict(
#         pad = 15,
#         thickness = 20,
#         line = dict(color = "black", width = 0.5),
#         label = sankey_data['labels'],
#         color = "#e0e0e0"
#         ),
#         link = dict(
#             source = sankey_data['sources'],
#             target = sankey_data['targets'],
#             value = sankey_data['values'],
#             color = sankey_data['colors']
#     ))])
#     fig.update_layout(
#         font=dict(size = 12, color = 'black'),
#         )

#     fig.write_image(outfile)





def draw_parallel_categories(resolution, outfile, metadata, show=True):
    targets = list(resolution.columns.values)

    last_group = {}
    for g, d in resolution.groupby(targets[-1]):
        last_group[g] = "<br>".join(d.index)

    resolution[targets[-1]] = resolution[targets[-1]].apply(lambda x: last_group[x])
    resolution['Categories'] = [metadata['values'].loc[i] for i in resolution.index]
    resolution = resolution.sort_values( ['Categories'] + targets[::-1])

    width = len(targets) * 150
    height = resolution.shape[0] * 35
    fig = px.parallel_categories(resolution, dimensions=targets, color="Categories", 
                    color_continuous_scale=px.colors.sequential.Plasma,
                    width=width, height=height
                    )
    # Make margins fit labels
    max_label_len = max([len(i) for i in resolution.index])

    fig.update_layout(coloraxis_showscale=False, margin=dict(l=20, r=max_label_len * 6, t=20, b=20))

    fig.write_image(outfile, engine="orca")
    if show:
        fig.show()