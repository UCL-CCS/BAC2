RESIDUE_DICTIONARY = {
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
    'DA': 'A',
    'DC': 'C',
    'DG': 'G',
    'DT': 'T'
}


def convert_aa3to1(aa3):
    try:
        aa1 = RESIDUE_DICTIONARY[aa3]
    except:
        aa1 = 'X'

    return aa1


def convert_aa1to3(aa1):
    map = {v: k for k, v in RESIDUE_DICTIONARY.items()}

    try:
        aa3 = map[aa1]
    except:
        aa3 = 'UNK'

    return aa3


def uniquify_list(lst):
    # Preserves order

    seen = set()
    seen_add = seen.add

    u_lst = [x for x in lst if not (x in seen or seen_add(x))]

    return u_lst