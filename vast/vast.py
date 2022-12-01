import pandas as pd
import os
import logging
import numpy as np
from pathlib import Path
import sys

path_root = Path(__file__)
sys.path.append(str(path_root))


from tree import get_snp_alignment
from patterns import (
    get_pattern, read_matrix_windows_and_get_patterns,
    calculate_scores)
from utils import (
    load_matrix, load_required_snps, pull_required_snps_from_matrix,
    get_final_snp_table,
    draw_resolution_ascii_graph
)


def optimization_loop(patterns, starting_pattern, delta_cutoff, max_targets):
    """
    Run optimization loop and return a list of chosen targets
    """
    logger = logging.getLogger('vast')
    max_targets = np.inf if max_targets == None else max_targets
    max_score = patterns.shape[1]**2 - patterns.shape[1]
    iteration = 0
    current_score = max_score
    delta = 0
    # Add constant value to starting pattern to easily add to pattern matrix
    # Constant needs to be order of magnitude higher than any value in pattern matrix
    const = 100 ** len(str(patterns.max()))
    result_pattern = starting_pattern
    go = True
    targets = []
    logger.info("Beginning optimiation loop.")
    while go:
        iteration += 1
        logger.info("Running iteration {}".format(iteration))
        # Add const to current pattern
        cur_pattern = np.multiply(result_pattern, const)
        # Add current pattern to all available patterns
        opt_pattern = np.add(cur_pattern, patterns)
        # Calculate entropy score for each pattern combined with current pattern
        scores = calculate_scores(opt_pattern)
        # TODO: Add more advanced filtering rather than just min score
        # 1. Get the best score that does not overlap with current target positions
        # 2. If there is a tie, pick a pattern that has higher resolution overall.
        # 3. Weigh by certain resolution targets (pass a target resolution) and run like a decision tree
        
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
                result_pattern,
                (max_score - current_score) / max_score))
        logger.debug("Current Pattern: {}".format(cur_pattern))
        logger.debug("Scores: {}".format(scores))
        logger.debug("Min Pattern: {}".format(min_pattern))
        logger.debug("Pattern Result: {}".format(result_pattern))
        logger.debug("Delta: {}".format(delta))
        targets.append(min_pattern)

        if iteration >= max_targets:
            go = False
            logger.info(
                "Ending optimization because max targets ({}) has been reached".format(max_targets)
            )
    return targets, ((max_score - current_score) / max_score) * 100



def run_vast(matrix, outfile, delta_cutoff, max_targets, window, offset, required_snps):
    logger = logging.getLogger('vast')
    logger.info("Loading SNP matrix: {}".format(os.path.abspath(matrix)))
    snp_matrix = load_matrix(matrix)
    logger.info("Found {0} genomes and {1} snps.".format(
        snp_matrix.shape[1], snp_matrix.shape[0]))
    if required_snps:
        # Pull required snps from matrix and set current pattern as starting pattern
        logger.info("Loading required SNPS")
        req_snps =  pull_required_snps_from_matrix(
            snp_matrix,
            load_required_snps(required_snps))
        logger.info("Found {} required SNPs".format(req_snps.shape[0]))
        starting_pattern = get_pattern(
            req_snps.values)
    else:
        starting_pattern = np.zeros(snp_matrix.shape[1], dtype=int)
    logger.debug("Starting Pattern: {}".format(starting_pattern))
    logger.info(
        "Finding patterns using a window size of {0} and an offset of {1}".format(
            window, offset
        ))
    patterns, snps = read_matrix_windows_and_get_patterns(snp_matrix, offset, window)
    logger.info("Found {0} patterns.".format(patterns.shape[0]))
    selected_patterns, score = optimization_loop(patterns, starting_pattern, delta_cutoff, max_targets)
    snp_results = get_final_snp_table(snps, selected_patterns, snp_matrix).reset_index()
    logger.info(
        "VaST identified {0} SNPs at {1} target regions that provide a resolution score of {2:.2f}%".format(
            snp_results.shape[0], len(selected_patterns), score
        ))
    if required_snps:
        req_snps['Target_ID'] = "REQ"
        req_snps = req_snps.reset_index().set_index(['Genome', 'Pos', 'Target_ID']).reset_index()
        snp_results = pd.concat([req_snps, snp_results])
    
    snp_results.to_csv(outfile, index=False, sep="\t")
    logger.info("SNP targets have been written to {}".format(outfile))
    

    
