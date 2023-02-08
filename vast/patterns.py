import logging
import pandas as pd
import numpy as np


try:
    from vast.utils import (
        pull_required_snps_from_matrix, load_required_snps)

except ImportError:
    from utils import (
        pull_required_snps_from_matrix, load_required_snps)

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
        snp_matrix.astype('U'), axis=1, return_inverse=True)[-1]



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
            start = position_list[pos]
            patterns.append(pattern_id)
            last_pos = position_list[pos]
    return patterns



def read_matrix_windows_and_get_patterns(matrix, offset, window):
    patterns = []
    snps = []
    for genome, d in matrix.groupby(level=0):
        d = d.droplevel([l for l in d.index.names if l != 'Pos'])
        d['window'] = windows(d.index.values, window, offset)
        # Drop any positions that are not valid (i.e. in offset regions)
        d = d[~(d['window'].isna())]
        for _, dd in d.groupby('window'):
            snps.append([(genome, p) for p in dd.index])
            patterns.append(get_pattern(dd.drop('window', axis=1).values))

    return np.stack(patterns), snps



def get_starting_pattern(matrix, required_snps):
    """
    Initialize resolution pattern either with required snps or
    at baseline (no resolution, all genomes are apart of the same group)
    """
    logger = logging.getLogger('vast')
    if required_snps:
        # Pull required snps from matrix and set current pattern as starting pattern
        logger.info("Loading required SNPS")
        required_snps =  pull_required_snps_from_matrix(
            matrix,
            load_required_snps(required_snps))
        logger.info("Found {} required SNPs".format(required_snps.shape[0]))
        starting_pattern = get_pattern(
            required_snps.values)
    else:
        starting_pattern = np.zeros(matrix.shape[1], dtype=int)
    logger.debug("Starting Pattern: {}".format(starting_pattern))
    return {
        'pattern': starting_pattern,
        'required_snps': required_snps}


def get_patterns(matrix, offset, window):
    """
    Search for snp patterns across matrix file in windows with a specified offset.
    Returns the differentiation pattern for each target and the SNPs (Genome, Loc)
    that are included in the target windows.
    """
    patterns, snps = read_matrix_windows_and_get_patterns(matrix, offset, window)
    return {
        'patterns': patterns,
        'snps': snps
    }



# def calculate_scores(opt_patterns):
#     scores = []
#     for row in opt_patterns:
#         _, counts = np.unique(row, return_counts=True)
#         scores.append(sum([c**2 - c for c in counts]))
#     return scores


def calculate_gini(opt_patterns, metadata):
    metadata=np.array(metadata)
    scores = []
    for row in opt_patterns:
        sizes = []
        gini_impurity = []
        for group in np.unique(row):
            group = np.where(row==group)[0]
            group_sz = len(group)
            sizes.append(group_sz)
            gini_impurity.append(
                1 - sum([(n/group_sz)**2 for n in
                         np.unique(metadata[group],
                                   return_counts=True)[1]]))        
        sizes = [s/len(row) for s in sizes]
        scores.append(sum([g*s for g, s in zip(gini_impurity, sizes)]))
    return scores

