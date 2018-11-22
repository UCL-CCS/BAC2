from parmed import residue

AMINO_RESIDUES = {
    'ALA': 'A',
    'ARG': 'R',
    'ASP': 'D',
    'ASN': 'N',
    'CYS': 'C',
    'GLU': 'E',
    'GLN': 'Q',
    'GLY': 'G',
    'HSD': 'H',
    'HIS': 'H',
    'HSE': 'H',
    'HSP': 'H',
    'ILE': 'I',
    'LEU': 'L',
    'LYS': 'K',
    'MET': 'M',
    'PHE': 'F',
    'PRO': 'P',
    'SER': 'S',
    'THR': 'T',
    'TRP': 'W',
    'TYR': 'Y',
    'VAL': 'V',
}

DNA_RESIDUES = {
    'DA': 'A',
    'DC': 'C',
    'DG': 'G',
    'DT': 'T'
}

RNA_RESIDUES = {
    'GUA': 'G',
    'CYT': 'C',
    'ADE': 'A',
    'THY': 'T',
    'URA': 'U',
    'A': 'A',
    'C': 'C',
    'G': 'G',
    'T': 'T',
    'U': 'U',
}

NUCLEIC_RESIDUES = {**RNA_RESIDUES, **DNA_RESIDUES}

RESIDUE_DICTIONARY = {**AMINO_RESIDUES, **NUCLEIC_RESIDUES}

REVERSE_REV_DICTIONARY = {v: k for k, v in RESIDUE_DICTIONARY.items()}

ANION_RESIDUES = list(residue.ANION_NAMES)

CATION_RESIDUES = list(residue.CATION_NAMES)

ION_RESIDUES = ANION_RESIDUES + CATION_RESIDUES

WATER_RESIDUES = residue.WATER_NAMES

POLYMER_BONDS = {
                 'protein': {'start': "C", 'end': "N"},
                 'rna': {'start': "O3'", 'end': "P"},
                 'dna': {'start': "O3'", 'end': "P"}
                }


def convert_resname3to1(aa3):

    if aa3 in RESIDUE_DICTIONARY:
        aa1 = RESIDUE_DICTIONARY[aa3]
    else:
        aa1 = 'X'

    return aa1


def convert_resname1to3(aa1):

    if aa1 in REVERSE_REV_DICTIONARY:
        aa3 = REVERSE_REV_DICTIONARY[aa1]
    else:
        aa3 = 'UNK'

    return aa3


def convert_resname_list(sequence, seq_format='letter'):
    
    if format not in ['fasta', 'letter']:
        raise RuntimeWarning('Resname list conversion to unknown format '
                             'requested - using "letter" as default')

    seq_res_letter = [convert_resname3to1(x) for x in sequence]

    if seq_format == 'letter':

        return seq_res_letter

    else:

        return ''.join(seq_res_letter)


def uniquify_list(lst):
    # Preserves order

    seen = set()
    seen_add = seen.add

    u_lst = [x for x in lst if not (x in seen or seen_add(x))]

    return u_lst
