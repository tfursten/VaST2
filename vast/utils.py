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
                    line = list(map(lambda x: x.upper(), line))
                    if len(set(line[1:]) - set(["A", "C", "T", "G"])):
                        raise ValueError("Encountered Invalid Base: Only A,C,T,G are valid in matrix.")
                    genome, position = split(line[0], "::")
                    line = [genome, position] + line[1:]
                    line = "\t".join(line)
                    outfile.write(line + "\n")


def replace_bases_with_numbers(matrix):
    matrix = matrix.apply(lambda x: x.str.upper())
    return matrix.replace({"A": 0, "C": 1, "G": 2, "T": 3})


def replace_numbers_with_bases(matrix):
    return matrix.replace({0: "A", 1: "C", 2:"G", 3:"T"})


def load_matrix(matrix):
    """
    Read in matrix and make sure to sort positions so we can calculate
    sliding windows. Replaces bases with numbers.
    """
    df = pd.read_csv(matrix, sep="\t", index_col=[0, 1]).sort_index()
    df.index.names = ['Genome', 'Pos']
    return replace_bases_with_numbers(df)


def load_required_snps(required_snps):
    return pd.read_csv(required_snps, sep="\t", header=None)

def pull_required_snps_from_matrix(matrix, required_snps):
    idx = [(row[0], row[1]) for _, row in required_snps.iterrows()]
    return matrix.loc[idx]

def get_final_snp_table(snps, selected_patterns, matrix):
    snp_df = pd.DataFrame([])
    for i in selected_patterns:
        required_snps = []
        for snp in snps[i]:
            required_snps.append(list(snp))
        d = pull_required_snps_from_matrix(
                matrix, pd.DataFrame(required_snps))
        cols = d.columns.values
        d['Target_ID'] = i
        snp_df = pd.concat([snp_df, d])

    return replace_numbers_with_bases(snp_df[['Target_ID'] + list(cols)])

def draw_resolution_ascii_graph(resolution, score):
    # get counts of group sizes
    groups, counts = np.unique(
        np.unique(resolution, return_counts=True)[1],
        return_counts=True)
    counts = np.array(np.divide(counts, counts.min()), dtype=int)
    chart = "\n" + "Resolution progress: {:.2f}%\n".format(score*100) + ("=" * 80) + "\n"
    chart += ("█" * int(np.floor(80*score))) + "\n"
    chart += ("=" * 80) + "\nGroup Size Distribution\n" 
    for group, count in zip(groups, counts):
        chart += ("█" * count) + " {0} group(s) of size {1}\n".format(count, group)
    chart += ("=" * 80)
    return chart
    