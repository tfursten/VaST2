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
        get_starting_pattern,
        calculate_gini)
    from vast.utils import (
        matrix_setup,
        get_final_snp_table,
        draw_resolution_ascii_graph, process_metadata
    )
except ImportError:
    from tree import (
        vast_target_tree,
        draw_parallel_categories)
    from patterns import (
        get_patterns,
        get_starting_pattern,
        calculate_gini)
    from utils import (
        matrix_setup,
        get_final_snp_table,
        draw_resolution_ascii_graph, process_metadata
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
    drop_duplicates, metadata, tree_outfile, figure_outfile):
    logger = logging.getLogger('vast')
    logger.info("Loading SNP matrix: {}".format(os.path.abspath(matrix)))
    snp_matrix = matrix_setup(matrix, exclude_snps, drop_duplicates)
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
    
    if tree_outfile:
        tree = vast_target_tree(snp_results, metadata)
        logger.info("Writing tree to file {}".format(tree_outfile))
        Phylo.write(tree, tree_outfile, "newick")
        logger.info(Phylo.draw_ascii(tree))
    if figure_outfile:
        draw_parallel_categories(results['resolution'], figure_outfile, metadata)


    
    

    
