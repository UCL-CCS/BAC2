import os
import sys
import json

from pathlib import Path

import freesasa

defaults_dirname = os.path.dirname(os.path.realpath(__file__))
default_config_filename = os.path.join(defaults_dirname, 'amber_config.txt')
default_atom_params_filename = os.path.join(defaults_dirname, 'wsas-params-wang2012.json')


class FreesasaRunner:
    """Wrapper to help run freesasa on a single PDB file

    Freesasa has a nice Python interface but some things don't work quite as
    needed for BAC, at least in Python 3. This wrapper is intended to handle
    these issues:
    1. File names need conversion to bytes when passed to freesasa
    2. By default should include HETATMS and hydrogens in analysis

    Parameters
    ----------

    config: str, optional
        Path to configuration file containing residue composition
        and atomic parameters - freesasa format.
    options: dict, optional
        Options to change how PDBs are parsed by freesasa.
    parameters: dict, optional
        Parameters to alter how freesasa computes surface area.

    Methods
    -------
    run(pdb)
        Run freesasa on input PDB file & return surface area results.

    """

    def __init__(self, **kwargs):
        """Wrapper for freesasa

        config: str, optional
            Path to configuration file containing residue composition
            and atomic parameters - freesasa format.
        options: dict, optional
            Options to change how PDBs are parsed by freesasa.
        parameters: dict, optional
            Parameters to alter how freesasa computes surface area.

        """

        default_options = {
            'hetatm': True,
            'hydrogen': True,
            # 'separate-chains' : False,
            'separate-models': True
        }

        default_parameters = {
            'algorithm': freesasa.LeeRichards,
            'probe-radius': freesasa.defaultParameters['probe-radius'],
            'n-points': freesasa.defaultParameters['n-points'],
            'n-slices': freesasa.defaultParameters['n-slices'],
            'n-threads': freesasa.defaultParameters['n-threads']
        }

        config_filename = kwargs.get('config', default_config_filename)
        config_filename = bytes(str(config_filename), 'utf-8')

        logfilename = kwargs.get('logfile', '/dev/null')

        # Hide warnings (as the load of multiple structures is two step and
        # extended config is not read in first step).
        freesasa.setVerbosity(1)

        self.classifier = freesasa.Classifier(config_filename)

        self.options = kwargs.get('options', default_options)

        self.parameters = kwargs.get('parameters', default_parameters)


    def run(self, pdb):
        """Run freesasa on provided PDB file

        Parameters
        ----------

        pdb: str
            Path to input PDB file

        Returns
        -------
        list
            SASA values for each atom of every model in the input PDB.

        """

        structure_array = freesasa.structureArray(bytes(pdb, 'utf-8'),
                                                  options=self.options,
                                                  classifier=self.classifier)

        results = []

        for s in structure_array:

            result = freesasa.calc(s)

            atom_areas = [result.atomArea(ndx) for ndx in range(s.nAtoms())]
            results.append(atom_areas)

        return results


def _create_freesasa_section_text(new_residues, sasa_atom_params):
    """
    Create text to add to freesasa configuration file to incorporate new residue.

    Parameters
    ----------
    new_residues : dict
        Non-standard residues to add to the freesasa config file.
        keys = residue names, values = atom name to type mapping (dict).
    sasa_atom_params: dict
        Maps atom type to properties needed by freesasa (radius and polarity).

    Returns
    -------
    atom_type_section : str
        Text to be added to freesasa config file atom type section.
    residue_section : str
        Text to be added to freesasa config file residue section.

    """

    atom_types = []

    # Create lines for residue section of format:
    # residue_name atom_name atom_type
    residue_section = ''

    for res_name, atom_to_type in new_residues.items():

        residue_section += '\n'

        for atom_name, atom_type in atom_to_type.items():
            residue_line = '{:s} {:s} {:s}\n'.format(res_name,
                                                     atom_name,
                                                     atom_type)

            atom_types.append(atom_type)

            residue_section += residue_line

    # Create lines for atom type section of format:
    # atom_type residue polarity
    atom_type_section = ''

    for atom_type in set(atom_types):

        atom_line = '{:s} {:.2f} {:s}\n'.format(atom_type,
                                                sasa_atom_params[atom_type]['radius'],
                                                sasa_atom_params[atom_type]['polarity'])

        atom_type_section += atom_line

    return atom_type_section, residue_section


def add_residues_freesasa_config_file(new_residues,
                                      new_filename,
                                      param_filename = default_atom_params_filename,
                                      orig_filename=default_config_filename):
    """
    Create a new freesasa config file that adds specified residue to the
    content of an existing copy.

    Parameters
    ----------
    new_residues : dict
        Non-standard residues to add to the freesasa config file.
        keys = residue names, values = atom name to type mapping (dict).
    new_filename: str
        Filename to be used for the updated freesasa config file.
    res_atom_to_type : dict
        Provides mapping from atom name to atom type.
    param_filename: str
        File containing atom type to properties (radius and polarity) mapping
        needed by freesasa.
    orig_filename: str
        Filename for the original freesasa config file.

    """

    with open(param_filename, 'r') as f:
        params = json.load(f)

    # Get text to add atom type and residue sections for the
    # residues being added to the config file
    (new_atom_types,
     new_residues) = _create_freesasa_section_text(new_residues,
                                                  params['params'])

    out_file = open(new_filename, 'w')

    with open(orig_filename) as input_config:

        for next_text in input_config:

            # Check for the start of non standard atom types
            # add new tyoes here
            if next_text.startswith('# extra'):
                next_text += new_atom_types

            out_file.write(next_text)

    # Insert new residue atom list
    out_file.write(new_residues)

    out_file.close()

    return

