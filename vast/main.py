import click
import logging
from pathlib import Path
import sys
import numpy as np

path_root = Path(__file__)
sys.path.append(str(path_root))

try:
    from vast.vast import run_vast, run_vast_resolution, run_vast_tree
    from vast.utils import (
        nasp_2_vast_format, run_get_snps_in_ranges,
        run_matrix_position_filter, run_write_snps_to_fasta)
except ImportError:
    from vast import run_vast, run_vast_resolution, run_vast_tree
    from utils import (
        nasp_2_vast_format, run_get_snps_in_ranges,
        run_matrix_position_filter, run_write_snps_to_fasta)

logging.basicConfig(format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
logger = logging.getLogger('vast')

@click.group()
@click.option('--debug/--no-debug', default=False)
def cli(debug):
    if debug:
        logger.setLevel('DEBUG')
        logger.debug("DEBUG LOGGING")
    else:
        logger.setLevel("INFO")


@cli.command(context_settings=dict(show_default=True))
@click.argument(
    'MATRIX',
    type=click.Path(exists=True, file_okay=True, dir_okay=False, allow_dash=False))
@click.argument(
    'OUTFILE', default="vast_matrix.tsv",
    type=click.Path(
        exists=False, file_okay=True, dir_okay=False, allow_dash=True, writable=True))
def nasp_format(matrix, outfile):
    """Convert a Nasp bestsnp tsv file into VaST format."""
    nasp_2_vast_format(matrix, outfile)

# TODO: allow other input format conversions like fasta, etc

@cli.command(context_settings=dict(show_default=True))
@click.argument(
    'MATRIX',
    type=click.Path(exists=True, file_okay=True, dir_okay=False, allow_dash=False))
@click.argument(
    'OUTFILE', default="vast_targets.tsv",
    type=click.Path(
        exists=False, file_okay=True, dir_okay=False, writable=True))
@click.option(
    '--metadata', '-c', 
    type=click.Path(exists=True, file_okay=True, dir_okay=False, allow_dash=False),
    help="File providing classification labels for each genome. Insead of maximizing resolution of all genomes (default), VaST will find targets that best differentiate between the categories. Two column tab separated file (no header, '#' for comments) should contain the strain/genome name and a category.")
@click.option(
    "--max-targets", '-m',
    type=click.IntRange(min=0), default=None,
    help="Stop search once MAX-TARGETS are selected. If DELTA value is provided, optimization will stop at whichever occurs first.")
@click.option(
    "--delta", '-d',
    default=0, type=click.FloatRange(min=0),
    help="Stop search when the resolution improvement is less than DELTA. If MAX_TARGETS is provided, optimization will stop at whichever occurs first.")
@click.option(
    "--window", '-w',
    default=50, type=click.IntRange(min=1),
    help="Size of scanning window for targets.")
@click.option(
    "--offset", '-o',
    default=1, type=click.IntRange(min=1),
    help="Set offset for target windows. Allows spacing between adjacent target regions.")
@click.option(
    "--required-snps", '-i',
    type=click.Path(exists=True, file_okay=True, dir_okay=False),
    help="File with list of SNPs that are required in the solution. Two column tab separated file should contain the reference genome id and the SNP positions (no header, '#' for comments). These will be labeled as 'REQ' in the output file.")
@click.option(
    "--exclude-snps", '-x',
    type=click.Path(exists=True, file_okay=True, dir_okay=False),
    help="File with a list of SNPs that should be excluded from the solution. Two column tab separated file should contain the reference genome id and the SNP positions (no header, '#' for comments)."
)
@click.option(
    "--drop-duplicates", '-r',
    default=False, is_flag=True,
    help="Speed up processing by removing duplicate columns (genomes) from SNP matrix (keep only first occurance)."
)
@click.option(
    "--fasta", '-f', default=None,
    type=click.Path(exists=False, file_okay=True, dir_okay=False, writable=True),
    help="Concatenate SNPs and write to fasta file."
)
@click.option(
    "--figure", '-s', default=None,
    type=click.Path(exists=False, file_okay=True, dir_okay=False, writable=True),
    help="Path to save Sankey diagram of resolution for chosen targets."
)
@click.option(
    "--resolution", '-z', default=None,
    type=click.Path(exists=False, file_okay=True, dir_okay=False, writable=True),
    help="Path to save the resolution pattern for each target (.tsv)."
)
def targets(
    matrix, outfile, delta, max_targets,
    window, offset, required_snps, exclude_snps,
    drop_duplicates, metadata, figure, fasta, resolution):
    """
    Run vast to find target regions to maximize strain resolution.\n
    VaST uses a greedy optimization algorithm to find a minimal number of target regions
    that will provide the maximal resolution of genomes. VaST takes a SNP matrix (columns
    are genomes and rows are SNP positions) and creates
    target regions by passing a sliding window across the SNP positions. Each collection of
    SNPs in a window has a different pattern in which it differentiates between the genomes.
    At each iteration, a pattern is chosen that
    provides the maximal increase in resolution of the genomes and it is added to the
    collection of targets.
    """
    logger.info("Running VaST")
    run_vast(
        matrix, outfile, delta, max_targets,
        window, offset, required_snps, exclude_snps,
        drop_duplicates, metadata, fasta, figure, resolution)


@cli.command(context_settings=dict(show_default=True))
@click.argument(
    'VAST_TARGET_MATRIX', 
    type=click.Path(exists=True, file_okay=True, dir_okay=False, allow_dash=False))
@click.option(
    '--metadata', '-c', 
    type=click.Path(exists=True, file_okay=True, dir_okay=False, allow_dash=False),
    help="File providing classification labels for each genome. Sankey diagram will be colored using these categories."
)
@click.option(
    "--figure", '-f', default=None,
    type=click.Path(exists=False, file_okay=True, dir_okay=False, writable=True),
    help="Path to save Sankey diagram of resolution for chosen targets."
)
@click.option(
    "--resolution", '-z', default=None,
    type=click.Path(exists=False, file_okay=True, dir_okay=False, writable=True),
    help="Path to save the resolution pattern for each target (.tsv)."
)
def resolution(vast_target_matrix, metadata, resolution, figure):
    """
    Calculate resolution from VaST target matrix (output from running vast targets).
    Output a Sankey diagram showing how each target splits up the collection of genomes
    and/or a table showing the same differentiation pattern (Diagram is not recommended
    when the number of genomes is high). 
    Input: A VaST target matrix with tab separated columns for 
    "Genome", "Pos", "Target_ID" followed
    by genome names. 
    """
    run_vast_resolution(vast_target_matrix, metadata, resolution, figure)


@cli.command(context_settings=dict(show_default=True))
@click.argument(
    'VAST_TARGET_MATRIX', 
    type=click.Path(exists=True, file_okay=True, dir_okay=False, allow_dash=False))
@click.option(
    '--metadata', '-c', 
    type=click.Path(exists=True, file_okay=True, dir_okay=False, allow_dash=False),
    help="File providing classification labels for each genome. Categories will be include in tree labels."
)
@click.option(
    "--tree", '-t', default=None,
    type=click.Path(exists=False, file_okay=True, dir_okay=False, writable=True),
    help="Write neighbor-joining tree using chosen targets to newick formated file."
)
def tree(vast_target_matrix, metadata, tree):              
    """
    Draw neighbor joining tree from VaST target matrix.  
    Input: A VaST target matrix with tab separated columns for 
    "Genome", "Pos", "Target_ID" followed
    by genome names. 
    """
    run_vast_tree(vast_target_matrix, metadata, tree)
#TODO: change tree drawing to use grapetree


@cli.command(context_settings=dict(show_default=True))
@click.argument(
    'SNP_MATRIX', 
    type=click.Path(exists=True, file_okay=True, dir_okay=False, allow_dash=False))
@click.argument(
    'RANGES',
    type=click.Path(exists=True, file_okay=True, dir_okay=False, allow_dash=False))
@click.argument(
    'OUTFILE',
    type=click.Path(exists=False, file_okay=True, dir_okay=False, writable=True))
@click.option(
    "--full-matrix", '-f',
    default=False, is_flag=True,
    help="By default, only SNP positions is output, using the full-matrix flag will return all SNP data in the range"
)
def get_snps_in_ranges(snp_matrix, ranges, outfile, full_matrix):
    """
    Given a RANGES tab separated file (with columns "Genome" for the reference
    genome name, "Start", and "End" indicating the start and end of the range),
    return a list of SNP positions from the provided SNP_MATRIX that fall within
    each range (inclusive). The output file is provided in a format that can be passed to the
    vast targets function as required or excluded snps.
    """
    run_get_snps_in_ranges(snp_matrix, ranges, outfile, full_matrix)

@cli.command(context_settings=dict(show_default=True))
@click.argument(
    'SNP_MATRIX', 
    type=click.Path(exists=True, file_okay=True, dir_okay=False, allow_dash=False))
@click.argument(
    'POSITIONS',
    type=click.Path(exists=True, file_okay=True, dir_okay=False, allow_dash=False))
@click.argument(
    'OUTFILE',
    type=click.Path(exists=False, file_okay=True, dir_okay=False, writable=True))
def matrix_position_filter(snp_matrix, positions, outfile):
    """
    Given a position file (with columns "Genome" for the reference
    genome name, and "Pos" for position of SNP in genome), pull only 
    the included positions from the provided SNP_MATRIX.
    """
    run_matrix_position_filter(snp_matrix, positions, outfile)


@cli.command(context_settings=dict(show_default=True))
@click.argument(
    'VAST_TARGET_MATRIX', 
    type=click.Path(exists=True, file_okay=True, dir_okay=False, allow_dash=False))
@click.argument(
    'OUTFILE',
    type=click.Path(exists=False, file_okay=True, dir_okay=False, writable=True))
def target_matrix_to_fasta(vast_target_matrix, outfile):
    """
    Convert a VaST Target Matrix (result from running vast target) into a fasta file
    by concatenating SNPs for each genome.
    """

    run_write_snps_to_fasta(vast_target_matrix, outfile)

if __name__ == '__main__':
    cli()



