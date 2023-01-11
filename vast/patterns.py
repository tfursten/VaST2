import numpy as np


def get_pattern(snp_matrix):
    """
    Given a matrix of snps for different genomes, return an 
    array that clusters genomes with the same genotype
    for the provided SNPs.
    Given: 
            [Genomes]
    [LOCI]  1  1  0  1  1
            2  2  2  1  1
            3  3  3  2  2
    Return:
    [0, 0, 1, 2, 2]
    """
    return np.unique(
        snp_matrix, axis=1, return_inverse=True)[-1]



def windows(position_list, window_size, offset):
    start = position_list[0]
    pattern_id = 0
    patterns = []
    last_pos = 0
    for pos in range(len(position_list)):
        if (position_list[pos] - start) < window_size:
            patterns.append(pattern_id)
            last_pos = position_list[pos]
        # position within offset
        elif (position_list[pos] - last_pos) < offset:
            patterns.append(None)
        else:
            pattern_id += 1
            start = pos
            patterns.append(pattern_id)
            last_pos = position_list[pos]
    return patterns



def read_matrix_windows_and_get_patterns(matrix, offset, window):
    patterns = []
    snps = []
    for genome, d in matrix.groupby(level=0):
        d = d.droplevel(0)
        d['window'] = windows(d.index.values, window, offset)
        # Drop any positions that are not valid (i.e. in offset regions)
        d = d[~(d['window'].isna())]
        for _, dd in d.groupby('window'):
            snps.append([(genome, p) for p in dd.index])
            patterns.append(get_pattern(dd.drop('window', axis=1).values))

    return np.stack(patterns), snps




# def read_matrix_windows_and_get_patterns(matrix, offset, window):
#     patterns = []
#     snps = []
#     current_genome = "NONE"
#     current_pos = 0
#     for i in range(matrix.shape[0]):
#         genome, pos = matrix.iloc[i].name
#         if (current_genome == genome) and ((pos - current_pos) < offset):
#             continue
#         else:
#             # TODO: Redo this, assume that matrix is sorted and just iterate
#             win = matrix.query(
#                 "Genome=='{genome}' and Pos >= {pos} and Pos < {window}".format(
#                     genome=genome, pos=pos, window=window + pos))
#             patterns.append(get_pattern(win.values))
#             snps.append(list(win.index.values))
#             # print(snps)
#             # print(list(win.index.values))
#             current_genome = genome
#             current_pos = pos
#     print(snps)
#     return np.stack(patterns), snps

def calculate_scores(opt_patterns):
    scores = []
    for row in opt_patterns:
        _, counts = np.unique(row, return_counts=True)
        scores.append(sum([c**2 - c for c in counts]))
    return scores