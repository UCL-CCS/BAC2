import re
from bac.builder.utils.sequence import convert_resname3to1
from bac.builder.utils.sequence import AMINO_RESIDUES, DNA_RESIDUES
from bac.builder.utils.sequence import RNA_RESIDUES, WATER_RESIDUES
from bac.builder.utils.sequence import ION_RESIDUES, POLYMER_BONDS


def get_residue_type(residue):
    """
    Return the rtype of the provided residue (based on the residue name).

    Parameters
    ----------
    residue : :class:`Residue <parmed.topologyobjects.Residue>`
        Residue on which to perform check.

    Returns
    -------
    str
        Residue type - 'protein', 'rna', 'dna', 'water' or 'unknown'
    """

    if residue.name in AMINO_RESIDUES.keys():
        residue_type = 'protein'
    elif residue.name in RNA_RESIDUES.keys():
        residue_type = 'rna'
    elif residue.name in DNA_RESIDUES.keys():
        residue_type = 'dna'
    elif residue.name in WATER_RESIDUES:
        residue_type = 'water'
    elif residue.name in ION_RESIDUES:
        residue_type = 'ion'
    else:
        residue_type = 'unknown'

    return residue_type


def check_residue_polymer_bonded(residue, bond_type='',
                                 directions=[1, -1], specific_partner=None):
    """
    Check if the supplier residue is recognisably bonded to known polymer
    type.

    Parameters
    ----------
    residue : :class:`Residue <parmed.topologyobjects.Residue>`
        Residue on which to perform check.
    bond_type : str
        Specify a particular polymer bond type to check (from those in
        `POLYMER_BONDS`)
    directions : list
        Directions along chain to check bonds, options only 1 or -1.
    specific_partner :class:`Residue <parmed.topologyobjects.Residue>` or None
        Another residue if the bond must be to a specified residue.

    Returns
    -------
    str
        Polymer type or None.

    """

    if not len([x for x in directions if x in [1, -1]]) == len(directions):
        raise RuntimeError("In check_residue_for_bond direction "
                           "keyword must be and 1 or -1.")

    atoms = residue.atoms

    if bond_type:
        check_bonds = {k: v for k, v in POLYMER_BONDS.items() if k == bond_type}
    else:
        check_bonds = POLYMER_BONDS

    for check_direction in directions:

        for polymer_type, bond_atoms in check_bonds.items():

            if check_direction == -1:
                bonded_atom_name = bond_atoms['end']
                partner_atom_name = bond_atoms['start']
            else:
                # check direction == 1 (known due to earlier check)
                bonded_atom_name = bond_atoms['start']
                partner_atom_name = bond_atoms['end']

            bonded_atom = [atom for atom in atoms if
                           atom.name == bonded_atom_name]

            if bonded_atom:

                bonded_atom = bonded_atom[0]

                partner_atom = [atom for atom in bonded_atom.bond_partners
                                if atom.name == partner_atom_name]

                if partner_atom:

                    partner_atom = partner_atom[0]

                    if specific_partner is not None:

                        if partner_atom.residue != specific_partner:
                            break

                    return polymer_type

    return None


def clean_residue_altlocs(struct, residue_idx,
                          altloc='A', blank_remaining=True):
    """
    Clean selected residue to leave only a single set of coordinates by
    removing all but the selected `altloc` and any common atoms.

    Parameters
    ----------
    struct : :class:`Structure <parmed.Structure>`
        Structure to check chain numbering in.
    residue_idx : int
        Index of the residue to clean.
    altloc : str
        Letter code of the altloc to retain.
    blank_remaining : bool
        Should we blank the altloc field on retained atoms.

    """

    keep = ['', altloc]

    mask = [1 if atom.altloc not in keep and atom.residue.idx == residue_idx
            else 0 for atom in struct]

    struct.strip(mask)

    if blank_remaining:
        for atom in struct.residues[residue_idx].atoms:
            atom.altloc = ''


def has_altloc(residue):
    """
    Check if residue contains altloc atoms.

    Parameters
    ----------
    residue : `Residue <parmed.topologyobjects.Residue>`
        Residue object containing to atoms to be checked.

    Returns
    -------
    bool
        Does the residue contain altlocs?

    """
    for atom in residue:
        if atom.altloc:
            return True
    return False


def get_altloc_residues(struct, idx_only=False):
    """
    Obtain list of all residues in `struct` containing altlocs.

    Parameters
    ----------
    struct : :class:`Structure <parmed.Structure>`
        Structure to check for altlocs.
    idx_only : bool
        Return index only (True) or residue objects (False = default).

    Returns
    -------
    list
        List containing either `parmed.topologyobjects.Residue` objects for all
        residues containing altlocs or their indices in the input structure.

    """

    have_altlocs = [residue for residue in struct.residues
                    if has_altloc(residue)]

    if not idx_only:
        return have_altlocs
    else:
        return [residue.idx for residue in have_altlocs]


def print_altloc_info(struct):
    """
    Print information on altlocs in `struct`

    Parameters
    ----------
    struct : :class:`Structure <parmed.Structure>`
        Structure about which to report altlocs.
    """

    print(f"idx\tname\tnumber\taltlocs")

    for residue in get_altloc_residues(struct):

        altlocs = [loc for loc in list(set([atom.altloc for atom in residue]))
                   if loc]

        print(f"{residue.idx}\t{residue.name}\t{residue.number}\t{altlocs}")


def residue_list_sequence(residues, seq_format='letter', gap_char='-'):
    """
    Get sequence information for selected chain in input structure.

    Parameters
    ----------
    residues : :class:`TrackedList <parmed.TrackedList>`
        Residue on which to perform check.
    seq_format : str
        Format for output - 'resname' (list of three letter reside codes),
        'letter' (list of single letter residue codes), 'fasta' (sting of
        single letter codes).
    gap_char : str
        Character to use in single letter format for gaps.

    Returns
    -------
    list or str
        Output of all residue codes in the format selected in seq_format.
    """

    if len(gap_char) > 1:
        raise RuntimeError("Must use a single character to represent "
                           "sequence gaps")

    if not residues:
        if seq_format == 'fasta':
            return ''
        else:
            return []

    if seq_format == 'resname':
        gap_code = gap_char*3
    else:
        gap_code = gap_char

    if seq_format == 'fasta':
        sequence = ''
    else:
        sequence = []

    last_chain = ''
    last_number = -1000

    for residue in residues:

        number_diff = residue.number - last_number
        if residue.chain == last_chain and number_diff > 1:

            for _ in range(number_diff - 1):
                if seq_format == 'fasta':
                    sequence += gap_code
                else:
                    sequence.append(gap_code)

        if seq_format in ['fasta', 'letter']:
            residue_letter = convert_resname3to1(residue.name)
            if seq_format == 'fasta':
                sequence += residue_letter
            else:
                sequence.append(residue_letter)
        else:
            sequence.append(residue.name)

        last_chain = residue.chain
        last_number = residue.number

    return sequence


def align_sequence_structure(residues, sequence):
    """
    Align sequence of residues in the selected chain of the structure with
    the provided sequence.

    Parameters
    ----------
    residues : :class:`TrackedList <parmed.TrackedList>`
        Residue on which to perform check.
    sequence : str
        FASTA-like sequence.

    Returns
    -------
    str
        Section of sequence that matches the structure.
    str
        Any residues in the input sequence before the matching section.
    str
        Any residues in the input sequence after the matching section.

    """

    alignment = ([], [], [])

    seq_regex = residue_list_sequence(residues, seq_format='fasta',
                                      gap_char='.')

    matches = re.search(seq_regex, sequence)

    if matches:
        pre_seq, post_seq = re.split(seq_regex, sequence)
        alignment = (matches[0], pre_seq, post_seq)

    return alignment
