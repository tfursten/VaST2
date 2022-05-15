from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment


def get_snp_alignment(matrix):
    """
    Given a SNP matrix with genome columns, create
    an alignment object.
    """
    seqs = []
    for c in matrix:
        seqs.append(SeqRecord(Seq("".join(matrix[c].values)), id=str(c)))

    return MultipleSeqAlignment(seqs)
    