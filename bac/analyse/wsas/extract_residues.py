import os
import parmed as pmd


def extract_parmtop_residue(filename):
    """
    Extract residue name and atom name/type mapping from input parmtop.
    Note: Only one residue must be present in the topology.

    Parameters
    ----------
    filename: str
        Filename of the input parmtop.

    Returns
    -------
    dict
        key = residue name, value = atom name to type mapping (dict).

    """

    res_top = pmd.load_file(str(filename))

    if len(res_top.parm_data['RESIDUE_LABEL']) == 1:

        res_name = res_top.parm_data['RESIDUE_LABEL'][0]
        res_atom_to_type = dict(zip(res_top.parm_data['ATOM_NAME'],
                                    res_top.parm_data['AMBER_ATOM_TYPE']))

    else:
        raise RuntimeError('Residue topology file must contain exactly one '
                           'residue - {:s}'.format(filename))

    return {res_name: res_atom_to_type}


def extract_prep_residue(filename):
    """
    Extract residue name and atom name/type mapping from input prep file.
    Note: Only one residue must be present.

    Parameters
    ----------
    filename: str
        Filename of the input prep file.

    Returns
    -------
    dict
        key = residue name, value = atom name to type mapping (dict).

    """

    res_name = None
    res_atom_to_type = {}

    with open(filename) as prep_file:

        line_no = 0

        for line in prep_file:

            if line_no == 4:
                res_name = line.split()[0]

            elif line_no > 9:

                if len(line.strip()) == 0:
                    break

                cols = line.split()

                res_atom_to_type[cols[1]] = cols[2]

            line_no += 1

    return {res_name: res_atom_to_type}


def extract_offlib_residue(filename):
    """
    Extract residue name(s) and atom name/type mapping from input prep file.

    Parameters
    ----------
    filename: str
        Filename of the input prep file.

    Returns
    -------
    dict
        keys = residue names, values = atom name to type mapping (dict).

    """

    lib_data = pmd.load_file(filename)

    residue_dict = {}

    for res_name, res_atoms in lib_data.items():
        residue_dict[res_name] = {atom.name: atom.type for atom in res_atoms}

    return residue_dict


def extract_mol2_residue(filename):
    """
    Extract residue name and atom name/type mapping from input prep file.
    Note: Only one residue must be present.

    Parameters
    ----------
    filename: str
        Filename of the input mol2 file.

    Returns
    -------
    dict
        keys = residue names, values = atom name to type mapping (dict).

    """

    mol2_data = pmd.load_file(filename)

    residue_dict = {mol2_data.name: {atom.name: atom.type for atom in mol2_data.atoms}}

    return residue_dict


def extract_residue(filename):
    """
    Extract residue name(s) and atom name/type mapping from an Amber format
    file containing residue template information. Filetype is guessed from
    extension.

    Parameters
    ----------
    filename: str
        Filename of the input residue template file.

    Returns
    -------
    dict
        keys = residue names, values = atom name to type mapping (dict).

    """

    extension = os.path.splitext(filename)[1]

    if extension in ['.off', '.lib']:
        residue = extract_offlib_residue(filename)

    elif extension in ['.top', '.prmtop']:

        residue = extract_parmtop_residue(filename)

    elif extension == '.prep':

        residue = extract_prep_residue(filename)

    elif extension == '.mol2':

        residue = extract_mol2_residue(filename)

    else:

        raise Exception('File {} is not a recognised Amber residue '
                        'template format'.format(filename))

    return residue


def extract_residues(filenames):
    """
    Extract residue name and atom name/type mapping from an Amber format
    file containing residue template information.

    Parameters
    ----------
    filenames: list
        Filenames of the input residue template files.

    Returns
    -------
    dict
        keys = residue names, values = atom name to type mapping (dict).

    """

    residues = {}

    for filename in filenames:
        residues.update(extract_residue(filename))

    return residues
