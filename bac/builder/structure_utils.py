import parmed as pmd
from bac.builder.utils.sequence import convert_resname_list


def chain_sequence(struct, chain_id, seq_format='letter'):
    """
    Get sequence information for selected chain in input structure.

    Parameters
    ----------
    struct : :class:`Structure <parmed.Structure>`
        Structure file from which to get chain sequence
    chain_id : str
        Selected chain
    seq_format: str
        Format for output - 'resname' (list of three letter reside codes),
        'letter' (list fo single letter residue codes), 'fasta' (sting of
        single letter codes).

    Returns
    -------
    list or str
        Output of all residue codes in the format selected in seq_format.
    """

    chain = struct.view[chain_id,:,:]

    sequence = [x.name for x in chain.residues]

    if seq_format == 'resname':
        return sequence

    else:

        return convert_resname_list(sequence, seq_format=seq_format)


def check_residue_type(struct, chain_id, resid, insertion):



    pass


def get_chain_gaps(struct, chain_id):

    pass


def reassign_chains():

    pass
