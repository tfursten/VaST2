import os
import logging
import pandas as pd
import numpy as np


def nasp_2_vast_format(matrix, outfile):
    last_row = 0
    split = str.split
    with open(matrix, 'r') as infile:
        with open(outfile, 'w') as outfile:
            for i, line in enumerate(infile):
                if i == 0:
                    header = split(line, "\t")
                    try:
                        last_row = header.index("#SNPcall")
                    except ValueError:
                        last_row = len(header)
                    header = header[1:last_row]
                    outfile.write("Genome\tPos\t" + "\t".join(header) + "\n")
                else:
                    line = split(line, "\t")
                    line = line[:last_row]
                    genome, position = split(line[0], "::")
                    line = [genome, position] + line[1:]
                    line = "\t".join(line)
                    outfile.write(line + "\n")


def load_matrix(matrix):
    """
    Read in SNP Matrix and make sure to sort positions so we can calculate
    sliding windows.
    """
    df = pd.read_csv(matrix, sep="\t", index_col=[0, 1], comment="#").sort_index()
    df.index.names = ['Genome', 'Pos']
    return df

def load_target_matrix(matrix):
    """
    Read in VaST Target matrix
    """
    df = pd.read_csv(matrix, sep="\t", index_col=[0,1,2], comment="#")
    return df


def load_required_snps(required_snps):
    return pd.read_csv(required_snps, sep="\t", header=None, comment="#", index_col=[0,1])

def pull_required_snps_from_matrix(matrix, required_snps):
    idx = required_snps.index
    idx_intersect = matrix.index.intersection(idx)
    return {
        'required_snps': matrix.loc[idx_intersect],
        'snp_matrix_without_required': matrix.drop(idx_intersect),
        'missing': idx.difference(idx_intersect)
        }

def filter_exluded_snps(snps, exclude_snps):
    """
    Remove exclude snps from snp matrix.
    """
    exclude_snps_df = pd.read_csv(
        exclude_snps, sep='\t', index_col=[0,1], header=None, comment="#")
    pre_filter = snps.shape[0]

    snps = snps.loc[~snps.index.isin(snps.index.intersection(exclude_snps.index))]
    post_filter = snps.shape[0]
    n_exclude = exclude_snps_df.shape[0]
    logging.info(
        "Filtered {0} SNPs out of {1} excluded SNPs".format(
        pre_filter - post_filter, n_exclude
    ))
    return snps


def get_final_snp_table(snps, selected_patterns, matrix, required_snps):
    if required_snps is not None:
        required_snps['Target_ID'] = "REQ"
        required_snps = required_snps.reset_index().set_index(
            ['Genome', 'Pos', 'Target_ID'])
    snp_dfs = []
    for n, i in enumerate(selected_patterns):
        selected_snps = []
        for snp in snps[i]:
            selected_snps.append(list(snp))
        d = pull_required_snps_from_matrix(
                matrix,
                pd.DataFrame(selected_snps).set_index([0,1])).get('required_snps')
        d['Target_ID'] = n
        d = d.set_index('Target_ID', append=True)
        d.index.set_names(['Genome', 'Pos', 'Target_ID'], inplace=True)
        snp_dfs.append(d)
    snp_results = pd.concat([required_snps] + snp_dfs)
    return snp_results

def draw_resolution_ascii_graph(resolution, score):
    # get counts of group sizes
    groups, counts = np.unique(
        np.unique(resolution, return_counts=True)[1],
        return_counts=True)
    norm_counts = np.array(
        np.interp(
            counts,
            (counts.min(), counts.max()),
            (counts.min(), min(counts.max(), 57))),
        dtype=int) # Normalize bar widths to fit within 80 characters.
    chart = "\n" + "Resolution progress: {:.2f}%\n".format(score*100) + ("=" * 80) + "\n"
    chart += ("█" * int(np.floor(80*score))) + "\n"
    chart += ("=" * 80) + "\nGroup Size Distribution\n" 
    for group, norm_count, count in zip(groups, norm_counts, counts):
        chart += ("█" * norm_count) + " {0} group(s) of size {1}\n".format(count, group)
    chart += ("=" * 80)
    return chart


def get_resolution(patterns, genomes):
    targets = ['Start'] + ["Target {}".format(i) for i in range(1, len(patterns))]
    return pd.DataFrame(
        patterns.T, index=genomes,
        columns = targets)


def process_metadata(matrix, metadata=None):
    """
    Get metadata categores for genomes.
    If no metadata is None, return an array with unique ids for each genome
    Otherwise, read metadata file and return classification ids for each
    genome in the same order as they are in the matrix.
    """
    genomes = matrix.columns
    if metadata is None:
        # If no metadata is provided treat each genome as a unique classification
        return {
            'values': pd.Series(np.arange(len(genomes)), index=genomes),
            'names': pd.Series(genomes)
        }
    metadata = pd.read_csv(metadata, sep='\t', index_col=0, header=None, comment="#")
    # ensure that the genomes in metadata match SNP matrix
    missing = np.setdiff1d(genomes, metadata.index)
    if len(missing):
        err = "Metadata missing for genomes: {}".format(", ".join(missing))
        raise(IndexError(err))
    # put metadata in same order as snp matrix
    metadata = metadata.loc[genomes]
    # Return unique ids for each category
    metadata_ids = np.unique(metadata[1], return_inverse=True)
    return {
        'values': pd.Series(metadata_ids[1], index=genomes),
        'names': pd.Series(metadata_ids[0])
    }


def get_metadata_from_genome(metadata, genome):
    return metadata['names'].loc[metadata['values'].loc[genome]]

def matrix_setup(matrix, exclude_snps, drop_duplicates):
    """
    Load SNP matrix from tsv file, and filter excluded SNPs
    """    
    snp_matrix = load_matrix(matrix)
    # drop duplicate genomes based on available SNPs
    # TODO: re-add duplicates to final matrix
    if drop_duplicates:
        snp_matrix = snp_matrix.T.drop_duplicates().T
    if exclude_snps:
        snp_matrix = filter_exluded_snps(snp_matrix, exclude_snps)
    return snp_matrix


def run_get_snps_in_ranges(snp_matrix, ranges, outfile, full_matrix=False):
    ranges = pd.read_csv(ranges, sep="\t", comment="#")
    ranges.columns = ['Genome', 'Start', 'End']
    useccols = None if full_matrix else [0,1] 
    snp_matrix = pd.read_csv(snp_matrix, sep="\t", usecols=useccols)
    snps_in_range = []
    for n, row in ranges.iterrows():
        snps_in_range.append(snp_matrix[
            (snp_matrix['Genome'] == row['Genome']) & 
            (snp_matrix['Pos'] >= row['Start']) &
            (snp_matrix['Pos'] <= row['End'])])
    snps_in_range = pd.concat(snps_in_range)
    snps_in_range.to_csv(outfile, sep="\t", header=False, index=False)


def run_matrix_position_filter(snp_matrix, positions, outfile):
    positions = load_required_snps(positions)
    snp_matrix = load_matrix(snp_matrix)
    filtered = pull_required_snps_from_matrix(snp_matrix, positions).get('required_snps')
    filtered.to_csv(outfile, sep="\t")


        

def write_snps_to_fasta(snp_matrix, fasta_outfile):
    with open(fasta_outfile, 'w') as out:
        for genome in snp_matrix:
            seq = "".join(snp_matrix[genome])
            out.write(
                ">{0}\n{1}\n".format(genome, seq))


def run_write_snps_to_fasta(vast_matrix, outfile):
    vast_matrix = load_target_matrix(vast_matrix)
    write_snps_to_fasta(vast_matrix, outfile)
