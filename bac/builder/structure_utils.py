import parmed as pmd
from bac.builder.utils.sequence import convert_resname_list


def chain_sequence(struct, chain_id, seq_format='letter'):

    chain = struct[chain_id,:,:]

    sequence = [x.name for x in chain.residues]

    if seq_format == 'resname':
        return sequence

    else:

        return convert_resname_list(sequence, seq_format=seq_format)


def check_residue_type(struct, chain_id, resid, insertion):

    pass


def get_chain_gaps(struct, chain_id):

    pass


