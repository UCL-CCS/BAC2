import numpy as np
import parmed as pmd
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
        Structure from which to get chain sequence
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

    chain = struct.view[chain_id, :, :]

    sequence = [x.name for x in chain.residues]

    if seq_format == 'resname':
        return sequence

    else:

        return convert_resname_list(sequence, seq_format=seq_format)


def scan_chain_type(struct, chain_id):
    """
    Lookup residue types (protein, water, etc.) for all residues in selected
    chain.

    Parameters
    ----------
    struct : :class:`Structure <parmed.Structure>`
        Structure to check residue types
    chain_id : str
        Selected chain

    Returns
    -------
    np.array
        Contains [residue index (in `struct`), residue type] pairs for each
        residue in selected chain.

    """
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


def check_residue_for_polymer_bond(struct, residue_idx,
                                   bond_type='protein', check_direction=1):
    """
    Determine if the selected residue in `struct` has a recognised polymer bond
    to the neighbouring residue.

    Parameters
    ----------
    struct : :class:`Structure <parmed.Structure>`
        Structure to check residue types
    residue_idx : int
        Index of the residue for which to check for a polymer bond
    bond_type : str
        Type of polymer bond to look for ('protein', 'dna', 'rna')
    check_direction : int
        Direction of bond to check along the polymer: either -1 toward lower
        residue numbers or 1 towards higher residue numbers.

    Returns
    -------
    bool
        True = residue is polymer bonded to neighbouring residue, False = no
        polymer bond to neighbouring residue
    """

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


def update_chain_type_assignment(struct, residue_type_assignment):
    """
    In the case where a chain contains multiple types we may have modified or
    artificial polymer residues incorrectly assigned. Here we check to see if
    at any change in assigned residue types the residue at the switch
    participate in polymer bonds and reassign type in
    `residue_type_assignment` accordingly.

    Parameters
    ----------
    struct : :class:`Structure <parmed.Structure>`
        Structure to check residue types
    residue_type_assignment : np.array
        Should contain [residue index (in `struct`), residue type] pairs for
        each residue in selected chain.

    Returns
    -------

    """

    # Find residues where type changes in  `residue_type_assignment`
    changes = np.where(residue_type_assignment[:-1, 1] != residue_type_assignment[1:, 1])[0]+1
    # We may check chain in either direction but should stop if find a second
    # type change or the end of the chain.
    search_ends = [0] + changes + residue_type_assignment[-1, 0]

    polymer_types = POLYMER_BONDS.keys()

    for idx in changes:

        if residue_type_assignment[idx, 1] != 'unknown':
            continue

        last_type = residue_type_assignment[idx + check_direction, 1]

        # If we are changing to a polymer type we should check backwards for
        # for residues that might change type, otherwise forwards.
        if last_type in polymer_types:
            check_direction = -1

        else:
            check_direction = 1

        while last_type in polymer_types:

            residue_idx = residue_type_assignment[idx, 0]

            polymer_bonded = check_residue_for_polymer_bond(struct,
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


def get_chain_gaps(struct, chain_id, return_idx=False):

    # TODO: Check this works

    chain = struct.view[chain_id, :, :]

    residue_numbers = np.array([[residue.number, residue.idx] for residue in chain])
    start_residue_idx = residue_numbers[np.where(np.diff(residue_numbers[:, 0]) != 1)][:, 1]

    if return_idx:

        return [(idx, idx + 1) for idx in start_residue_idx]

    else:

        residues = struct.residues

        return [(residues[idx].number, residues[idx + 1].number) for idx in start_residue_idx]


def clean_residue_altlocs(structure, residue_idx, altloc='A', blank_remaining=True):
    keep = ['', altloc]

    mask = [1 if atom.altloc not in keep and atom.residue.idx == residue_idx else 0 for atom in structure]

    structure.strip(mask)

    if blank_remaining:
        for atom in structure.residues[residue_idx].atoms:
            atom.altloc = ''


def has_altloc(residue):
    for atom in residue:
        if atom.altloc:
            return True
    return False


def get_altloc_residues(structure, idx_only=False):
    have_altlocs = []

    for residue in structure.residues:
        if has_altloc(residue):
            if idx_only:
                have_altlocs.append(residue.idx)
            else:
                have_altlocs.append(residue)

    return have_altlocs


def print_altloc_info(structure):

    print(f"idx\tname\tnumber\taltlocs")

    for residue in get_altloc_residues(structure):

        altlocs = [loc for loc in list(set([atom.altloc for atom in residue])) if loc]

        print(f"{residue.idx}\t{residue.name}\t{residue.number}\t{altlocs}")


def reassign_chains():

    pass