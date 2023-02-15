import pandas as pd
import os
import logging
import numpy as np
from Bio import Phylo
from pathlib import Path
import sys
path_root = Path(__file__)
sys.path.append(str(path_root))

try:
    from vast.tree import (
        vast_target_tree,
        draw_parallel_categories)
    from vast.patterns import (
        get_patterns,
        get_pattern,
        get_starting_pattern,
        calculate_gini)
    from vast.utils import (
        matrix_setup,
        get_final_snp_table,
        get_resolution,
        draw_resolution_ascii_graph,
        process_metadata
    )
except ImportError:
    from tree import (
        vast_target_tree,
        draw_parallel_categories)
    from patterns import (
        get_patterns,
        get_pattern,
        get_starting_pattern,
        calculate_gini)
    from utils import (
        matrix_setup,
        get_final_snp_table,
        get_resolution,
        draw_resolution_ascii_graph,
        process_metadata
    )



def optimization_loop(
    patterns, starting_pattern, delta_cutoff, max_targets, metadata):
    """
    Run optimization loop and return a list of chosen targets
    """
    logger = logging.getLogger('vast')
    max_targets = np.inf if max_targets == None else max_targets
    iteration = 0
    current_score = 1
    delta = 0
    # Add constant value to starting pattern to easily add to pattern matrix
    # Constant needs to be order of magnitude higher than any value in pattern matrix
    const = 100 ** len(str(patterns.max()))
    result_pattern = starting_pattern
    go = True
    targets = []
    resolution = []
    logger.info("Beginning optimiation loop.")
    while go:
        iteration += 1
        logger.info("Running iteration {}".format(iteration))
        # Add const to current pattern
        cur_pattern = np.multiply(result_pattern, const)
        # Add current pattern to all available patterns
        opt_pattern = np.add(cur_pattern, patterns)
        # Calculate entropy score for each pattern combined with current pattern
        scores = calculate_gini(opt_pattern, metadata)      
        # Find minimum score
        min_score = min(scores)
        # Find index of pattern with min score
        min_pattern = np.argmin(scores)
        # Update delta 
        delta = (current_score - min_score)
        # Stop if score does not improve by more than the cutoff
        if delta <= delta_cutoff:
            logger.info(
                "Ending optimization because score improvement is less than cutoff: {0} <= {1}".format(
                    delta, delta_cutoff))
            break
        current_score = min_score
        # update result pattern
        result_pattern = np.unique(opt_pattern[min_pattern], return_inverse=True)[-1]
        logger.info(
            draw_resolution_ascii_graph(
                result_pattern, 1 - current_score))
        logger.debug("Current Pattern: {}".format(cur_pattern))
        logger.debug("Scores: {}".format(scores))
        logger.debug("Min Pattern: {}".format(min_pattern))
        logger.debug("Pattern Result: {}".format(result_pattern))
        logger.debug("Delta: {}".format(delta))
        targets.append(min_pattern)
        resolution.append(list(result_pattern))


        if iteration >= max_targets:
            go = False
            logger.info(
                "Ending optimization because max targets ({}) has been reached".format(max_targets)
            )
    return {
        'targets': targets,
        'score': (1 - current_score) * 100,
        'resolution': np.array(resolution)
    }


def run_vast(
    matrix, outfile, delta_cutoff, max_targets,
    window, offset, required_snps, exclude_snps,
    drop_duplicates, metadata, tree_outfile,
    figure_outfile, resolution_outfile):
    logger = logging.getLogger('vast')
    logger.info("Loading SNP matrix: {}".format(os.path.abspath(matrix)))
    snp_matrix = matrix_setup(matrix, exclude_snps, drop_duplicates)
    genomes = snp_matrix.columns.values
    logger.info("Found {0} genomes and {1} snps.".format(
        snp_matrix.shape[1], snp_matrix.shape[0]))
    starting_pattern = get_starting_pattern(snp_matrix, required_snps)
    logger.info(
        "Finding patterns using a window size of {0} and an offset of {1}".format(
            window, offset
    ))
    patterns = get_patterns(snp_matrix, offset, window)
    logger.info("Found {0} patterns.".format(patterns['patterns'].shape[0]))
    try:
        metadata = process_metadata(snp_matrix, metadata)
    except IndexError as err:
        logger.error(err)
    results = optimization_loop(
        patterns['patterns'], starting_pattern['pattern'],
        delta_cutoff, max_targets, metadata['values'])
    snp_results = get_final_snp_table(
        patterns['snps'], results['targets'], snp_matrix,
        starting_pattern['required_snps']).reset_index()
    n_snps = snp_results[~ snp_results['Target_ID'].isin(["REQ"])].shape[0]
    logger.info(
        "VaST identified {0} SNPs at {1} target regions that provide a resolution score of {2:.2f}%".format(
            n_snps, len(results['targets']), results['score']
        ))

    snp_results.to_csv(outfile, index=False, sep="\t")
    logger.info("SNP targets have been written to {}".format(outfile))
    # Get resolution combining starting pattern and selected target patterns
    resolution = get_resolution(
        np.vstack((starting_pattern['pattern'], results['resolution'])), genomes)
    if resolution_outfile:
        logger.info("Writing resolution table to file {}".format(resolution_outfile))
        resolution.to_csv(resolution_outfile, sep="\t")
    if tree_outfile:
        tree = vast_target_tree(
            snp_results.set_index(
                ['Genome', 'Pos', 'Target_ID']),
            metadata)
        logger.info("Writing tree to file {}".format(tree_outfile))
        Phylo.write(tree, tree_outfile, "newick")
        logger.info(Phylo.draw_ascii(tree))
    if figure_outfile:
        logger.info("Writing resolution file to file {}".format(figure_outfile))
        draw_parallel_categories(resolution, figure_outfile, metadata)
    
    

def run_vast_resolution(vast_targets, metadata, resolution_outfile, figure_outfile):
    logger = logging.getLogger('vast')
    logger.info("Loading VaST Target Matrix: {}".format(os.path.abspath(vast_targets)))
    target_df = pd.read_csv(
        vast_targets, sep="\t", comment="#",
        index_col=['Genome', 'Pos', 'Target_ID'], converters={'Target_ID': lambda x: -1 if x=="REQ" else x})
    logger.info("Found {0} genomes, {1} targets, and {2} snps.".format(
        target_df.shape[1], target_df.groupby(level=2).ngroups , target_df.shape[0]))    
    try:
        metadata = process_metadata(target_df, metadata)
    except IndexError as err:
        logger.error(err)

    const = 10 ** len(str(target_df.shape[1]))
    patterns = []
    current_pattern = None
    for n, d in target_df.groupby(level=2):
        if current_pattern is None:
            if n != -1:
                # if there are no required snps, start with
                # no resolution.
                patterns.append(np.zeros(target_df.shape[1], dtype=int))
            current_pattern = get_pattern(d.values)
            patterns.append(np.array(current_pattern))
        else:
            prev_pattern = np.multiply(current_pattern, const)
            # Add current pattern to all available patterns
            dd = get_pattern(d.values)
            next_pattern = np.unique(
                np.add(prev_pattern, get_pattern(d.values)),
                return_inverse=True)[1]
            current_pattern = next_pattern
            patterns.append(current_pattern)
    patterns = np.vstack(patterns)
    resolution = get_resolution(patterns, list(target_df.columns.values))
    if resolution_outfile:
        logger.info("Writing resolution table to file {}".format(resolution_outfile))
        resolution.to_csv(resolution_outfile, sep="\t")

    if figure_outfile:
        logger.info("Writing resolution file to file {}".format(figure_outfile))
        draw_parallel_categories(resolution, figure_outfile, metadata)
    

def run_vast_tree(vast_targets, metadata, tree_outfile):
    logger = logging.getLogger('vast')
    logger.info("Loading VaST Target Matrix: {}".format(os.path.abspath(vast_targets)))
    target_df = pd.read_csv(
        vast_targets, sep="\t", comment="#",
        index_col=['Genome', 'Pos', 'Target_ID'], converters={'Target_ID': lambda x: -1 if x=="REQ" else x})
    logger.info("Found {0} genomes, {1} targets, and {2} snps.".format(
        target_df.shape[1], target_df.groupby(level=2).ngroups , target_df.shape[0]))    
    try:
        metadata = process_metadata(target_df, metadata)
    except IndexError as err:
        logger.error(err)

    if tree_outfile:
        tree = vast_target_tree(target_df, metadata)
        logger.info("Writing tree to file {}".format(tree_outfile))
        Phylo.write(tree, tree_outfile, "newick")
    logger.info(Phylo.draw_ascii(tree))
