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



def read_matrix_windows_and_get_patterns(matrix, offset, window, start=0):
    patterns = []
    snps = []
    current_genome = "NONE"
    current_pos = 0
    current_pattern = []
    current_snps = []
    for i in range(start, matrix.shape[0]):
        row = matrix.iloc[i]
        genome, pos = row.name
        if current_genome == "NONE":
            current_genome = genome
            current_pos = pos
        if (pos - current_pos < window) and (current_genome == genome):
            current_pattern.append(row.values)
            current_snps.append(row.name)
            current_pos = pos
        else:
            # add previous patterns
            if len(current_pattern):
                patterns.append(get_pattern(current_pattern))
                snps.append(current_snps)
                current_pattern = []
                current_snps = []
            # skip if within the offset from previous pattern
            if pos - current_pos < offset:
                continue
            else:
                # Add pattern and update positions
                current_pattern = [row.values]
                current_snps = [row.name]
                current_genome = genome
                current_pos = pos
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