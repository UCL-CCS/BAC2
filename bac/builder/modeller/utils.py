import textwrap


def write_ali_file(filename, target_seq, structure_seq,
                   structure_id='pdb', start_resid=1, end_resid=-1, chain='?'):
    """
    Write a single chain Modeller (PIR format) alignment file.

    Parameters
    ----------
    filename : str
        Filename for output
    target_seq : str
        Fasta style list of single letter residue codes of the sequence to
        be modelled.
    structure_seq : str
        Fasta style list of single letter residue codes of the sequence in the
        existing structure (with gaps marked by -s)
    structure_id : str
        A string to identify the input structure
    start_resid : str
        First residue number in the chain being modelled
    end_resid : str
        End residue number of the chain being modelled
    chain : str
        Chain being modelled

    """
    with open(filename, 'w') as output_file:

        output_file.write(f'>P1;{structure_id}\n')
        output_file.write(f'structureX:{structure_id}:{start_resid}:{chain}:{end_resid}:{chain}::::\n')
        for line in textwrap.wrap(structure_seq + '*', width=75, break_on_hyphens=False):
            output_file.write(line + '\n')
        output_file.write(f'>P1;{structure_id}_fill\n')
        output_file.write('sequence:::::::::\n')
        for line in textwrap.wrap(target_seq + '*', width=75, break_on_hyphens=False):
            output_file.write(line + '\n')