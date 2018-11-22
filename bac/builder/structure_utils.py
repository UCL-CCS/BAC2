import parmed as pmd
import numpy as np
from bac.builder.utils.sequence import convert_resname_list
from bac.builder.utils.sequence import AMINO_RESIDUES, DNA_RESIDUES
from bac.builder.utils.sequence import RNA_RESIDUES, WATER_RESIDUES
from bac.builder.utils.sequence import ION_RESIDUES, POLYMER_BONDS


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


def scan_chain(struct, chain_id):
    residue_type = []
    chain = struct.view[chain_id, :, :]

    for res in chain.residues:

        if res.name in AMINO_RESIDUES.keys():
            residue_type.append([res.idx, 'protein'])
        elif res.name in RNA_RESIDUES.keys():
            residue_type.append([res.idx, 'rna'])
        elif res.name in DNA_RESIDUES.keys():
            residue_type.append([res.idx, 'dna'])
        elif res.name in WATER_RESIDUES:
            residue_type.append([res.idx, 'water'])
        elif res.name in ION_RESIDUES:
            residue_type.append([res.idx, 'ion'])
        else:
            residue_type.append([res.idx, 'unknown'])

    return np.array(residue_type, dtype=object)


def check_residue_for_bond(struct, residue_idx,
                           bond_type='protein', check_direction=1):

    if bond_type not in POLYMER_BONDS:
        raise RuntimeError(f"Invalid polymer bond type, "
                           f"must be one of {polymer_bonds.keys()}")

    if check_direction == -1:
        bonded_atom_name = POLYMER_BONDS[bond_type]['end']
        partner_atom_name = POLYMER_BONDS[bond_type]['start']
    elif check_direction == 1:
        bonded_atom_name = POLYMER_BONDS[bond_type]['start']
        partner_atom_name = POLYMER_BONDS[bond_type]['end']
    else:
        raise RuntimeError("In check_residue_for_bond direction "
                           "keyword must be and 1 or -1.")

    res = struct.residues[residue_idx]

    bond_atom_idx = [x.idx for x in res.atoms if x.name == bonded_atom_name]

    if bond_atom_idx:

        bond_atom_idx = bond_atom_idx[0]

        bond_atom = struct.atoms[bond_atom_idx]

        bonded_atom = [x.idx for x in bond_atom.bond_partners
                       if x.name == partner_atom_name]

        if bonded_atom:
            return True

    return False


def check_chain_type_assignment(struct, residue_type_assignment):

    # Search `residue_type_assignment` to find residues where type changes
    changes = np.where(residue_type_assignment[:-1, 1] != residue_type_assignment[1:, 1])[0]+1
    search_ends = [0] + changes + residue_type_assignment[-1, 0]

    polymer_types = POLYMER_BONDS.keys()

    for idx in changes:

        if last_type in polymer_types:
            check_direction = -1

        else:
            check_direction = 1

        last_type = residue_type_assignment[idx + check_direction, 1]

        while last_type in polymer_types:

            residue_idx = residue_type_assignment[idx,0]

            polymer_bonded = check_residue_for_bond(struct,
                                                    residue_idx,
                                                    bond_type=last_type,
                                                    check_direction=check_direction)

            if polymer_bonded:

                residue_type_assignment[idx, 1] = last_type

                idx -= check_direction

                if idx in search_ends:
                    break

            else:
                break


def get_chain_gaps(struct, chain_id):

    pass


def reassign_chains():

    pass
