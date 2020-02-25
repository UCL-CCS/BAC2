import os

import freesasa

from .extract_residues import extract_residue

# Defaults

_DEFAULT_OPTIONS = {
    'hetatm': True,
    'hydrogen': True,
    # 'separate-chains' : False,
    'separate-models': True
}

_DEFAULT_PARAMETERS = {
    'algorithm': freesasa.LeeRichards,
    # 'probe-radius': freesasa.defaultParameters['probe-radius'],
    # 'n-points': freesasa.defaultParameters['n-points'],
    # 'n-slices': freesasa.defaultParameters['n-slices'],
    # 'n-threads': freesasa.defaultParameters['n-threads']
}


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

    def __init__(self, config, wsas_params, tmp_dir, nonstandard_residue_files, ligand_topology, options=None, parameters=None):
        """Wrapper for freesasa

        config: str
            Path to configuration file containing residue composition
            and atomic parameters - freesasa format.
        options: dict, optional
            Options to change how PDBs are parsed by freesasa.
        parameters: dict, optional
            Parameters to alter how freesasa computes surface area.

        """

        # Hide warnings (as the load of multiple structures is two step and
        # extended config is not read in first step).
        freesasa.setVerbosity(1)

        config = self._update_sasa_config(config, wsas_params, tmp_dir, nonstandard_residue_files, ligand_topology)

        self.classifier = freesasa.Classifier(bytes(str(config), 'utf-8'))

        self.options = options or _DEFAULT_OPTIONS

        self.parameters = parameters or _DEFAULT_PARAMETERS

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

        structure_array = freesasa.structureArray(bytes(pdb, 'utf-8'), options=self.options, classifier=self.classifier)

        results = []

        for s in structure_array:

            result = freesasa.calc(s)

            atom_areas = [result.atomArea(ndx) for ndx in range(s.nAtoms())]
            results.append(atom_areas)

        return results

    def _update_sasa_config(self, config, parameters, tmp_dir, nonstandard_residue_files, ligand_topology):
        """
        Add non-standard residues (including the ligand if a topology is
        provided for it) to the freesasa config file.

        Parameters
        ----------

        Notes
        -----
        Edited config files is saved in self.tmp_dir and
        self.freesasa_config_file is updated to reflect this.

        Returns
        -------

        """
        files_to_add = nonstandard_residue_files

        if ligand_topology:
            files_to_add.append(ligand_topology)

        residues_to_add = {}

        for filename in files_to_add:
            residues_to_add.update(extract_residue(filename))

        if residues_to_add:

            sasa_config = os.path.join(tmp_dir, 'system_sasa.config')

            self._add_residues_freesasa_config_file(residues_to_add, sasa_config, parameters, orig_filename=config)

            return sasa_config

        return config

    @staticmethod
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

    def _add_residues_freesasa_config_file(self, new_residues, new_filename, atom_params, orig_filename):
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
        atom_params: dict
            Radius and polarity information for each atom type.
        orig_filename: str
            Filename for the original freesasa config file.

        """

        # Get text to add atom type and residue sections for the
        # residues being added to the config file
        (new_atom_types, new_residues) = self._create_freesasa_section_text(new_residues, atom_params)

        with open(new_filename, 'w') as out_file, open(orig_filename) as input_config:
            [out_file.write(l+new_atom_types if l.startswith('# extra') else l) for l in input_config]
            out_file.write(new_residues)
