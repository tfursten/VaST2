import click
import logging
from pathlib import Path
import sys

path_root = Path(__file__)
sys.path.append(str(path_root))

from vast import run_vast
from utils import nasp_2_vast_format
import numpy as np

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


@cli.command()
@click.argument(
    'MATRIX',
    type=click.Path(exists=True, file_okay=True, dir_okay=True, allow_dash=True))
@click.argument(
    'OUTFILE', default="vast_matrix.tsv",
    type=click.Path(
        exists=False, file_okay=True, dir_okay=False, allow_dash=True, writable=True))
def nasp_format(matrix, outfile):
    """Convert a Nasp bestsnp tsv file into VaST format."""
    nasp_2_vast_format(matrix, outfile)

@cli.command()
@click.argument(
    'MATRIX',
    type=click.Path(exists=True, file_okay=True, dir_okay=True, allow_dash=True))
@click.argument(
    'OUTFILE', default="vast_targets.tsv",
    type=click.Path(
        exists=False, file_okay=True, dir_okay=False, writable=True))
@click.option(
    "--max-targets", '-m',
    type=click.IntRange(min=0),
    help="Stop search once MAX-TARGETS are selected. If DELTA value is provided, optimization will stop at whichever occurs first.")
@click.option(
    "--delta", '-d',
    default=0, type=click.FloatRange(min=0),
    help="Stop search when the resolution improvement is less than DELTA. If MAX_TARGETS is provided, optimization will stop at whichever occurs first."
)
@click.option(
    "--window", '-w',
    default=50, type=click.IntRange(min=0),
    help="Size of scanning window for targets.")
@click.option(
    "--offset", '-o',
    default=1, type=click.IntRange(min=1),
    help="Set offset for window overlap. An offset of 1 will check all possible windows but this is less efficient with many SNPs and can result in overlapping target regions.")
@click.option(
    "--required-snps", '-r',
    type=click.Path(exists=True, file_okay=True, dir_okay=False),
    help="File with list of SNPs that are required in the solution. Two column tab separated file should contain the genome id and the SNP positions")
def targets(matrix, outfile, delta, max_targets, window, offset, required_snps):
    """
    Run vast to find target regions.\n
    VaST uses a greedy optimization algorithm to find a minimal number of target regions
    that will provide the maximal resolution of genomes. VaST takes a SNP matrix (columns
    are genomes and rows are SNP genotypes) and creates
    target regions by passing a sliding window across the SNP positions. The collection
    of SNPs in each window are used to create a pattern for that locus which identifies
    which genomes have identical genotypes. At each iteration, a pattern is chosen that
    provides the maximal increase in resolution of the genomes and it is added to the
    collection of targets.
    """
    logger.info("Running VaST")
    run_vast(matrix, outfile, delta, max_targets, window, offset, required_snps)
    



# Add a number of filter functions to remove genomes, regions etc.
# Add a function to find a similar pattern as a given range (if we can't design primers for a certain target)
# Add viz
# Make an easy tool to convert genome ranges to a list of required positions.
# 



                





if __name__ == '__main__':
    cli()



