import os
import json
from pathlib import Path
from typing import List
from enum import Enum
import tempfile
import numpy as np

import parmed as pmd
import mdtraj as md
import pandas as pd

from bac.utils.decorators import advanced_property, pathlike
from .extract_residues import extract_residue
from .freesasa_utils import FreesasaRunner, add_residues_freesasa_config_file, default_config_filename

defaults_dirname = os.path.dirname(os.path.realpath(__file__))
default_params_filename = os.path.join(defaults_dirname, 'wsas-params-wang2012.json')

class Component(Enum):
    complex = 'complex'
    ligand = 'ligand'
    receptor = 'receptor'


class Wsas:

    def __init__(self, **kwargs):

        self.temperature = kwargs.get('temperature', 300)
        self.slice = slice(kwargs.get('first_frame', 0), kwargs.get('last_frame', -1), kwargs.get('stride', 1))
        self.ligand_filter = kwargs.get('ligand_filter', None)
        self.nonstandard_residue_files = kwargs.get('nonstandard_residue_files', [])

        default_solvent = [ "WAT", "HOH", "'Cl.*'", "CIO", "'Cs+'", "IB", "'K.*'",
                            "'Li+'", "'MG.*'", "'Na+'", "'Rb+'", "CS", "RB", "NA",
                            "F",  "ZN"]
        self.solvent_residues = kwargs.get('solvent_residues', default_solvent)

        self.topology = kwargs.get('topology')
        self.ligand_topology = kwargs.get('ligand_topology', None)
        self.trajectories: List[Path] = kwargs.get('trajectories', [])
        self.component = kwargs.get('component')
        self.tmp_dir = tempfile.mkdtemp()

        self.areas = {}
        self.wsas = {}
        self.energies = {}

        parameter_file = kwargs.get('parameter_file', default_params_filename)

        with open(parameter_file, 'r') as f:
            self.parameters = json.load(f)

        self.freesasa_config_file = kwargs.get('config_file', default_config_filename)

        if not self.ligand_filter and self.ligand_topology:

            lig_res = extract_residue(self.ligand_topology)

            self.ligand_filter = 'resname {:s}'.format(list(lig_res.keys())[0])


    @advanced_property(type=str)
    def topology(self): pass

    @advanced_property(type=Component, default=Component.complex)
    def component(self): pass

    @pathlike(default='.')
    def output(self): pass

    def calc_surface_areas(self, energy=True):

        self.output.mkdir(exist_ok=True, parents=True)

        topology_filename = self.topology
        top = pmd.load_file(topology_filename)

        self.update_sasa_config()

        sasa_config = self.freesasa_config_file
        sasa_calculator = FreesasaRunner(config=sasa_config)
        atom_selections = None

        results = {}

        for trajectory in self.trajectories:

            traj =  md.load(trajectory, top=topology_filename)

            # If this is the first trajectory setup dataframes to store data
            if atom_selections is None:

                atom_selections = self.create_component_selections(traj)

                for component, atom_list in atom_selections.items():
                    results[component] = pd.DataFrame()
                    results[component]['resid'] = [top.atoms[x].residue.number for x in atom_list]
                    results[component]['residue_name'] = [top.atoms[x].residue.name for x in atom_list]
                    results[component]['atom_name'] = [top.atoms[x].name for x in atom_list]
                    results[component]['atom_type'] = [top.atoms[x].type for x in atom_list]

            for component, atom_list in atom_selections.items():

                trajectory_pdb = '{}/{}-traj.pdb'.format(self.tmp_dir, component)
                comp_traj = traj[self.slice].atom_slice(atom_list)
                comp_traj.save_pdb(trajectory_pdb)

                # Compute atomic SASA
                atom_areas = pd.DataFrame(sasa_calculator.run(trajectory_pdb))

                # SASA not prone to large changes - so use average
                # TODO: option to get per frame values as output
                avg_areas = atom_areas.mean()
                results[component][trajectory] = avg_areas

        self.areas = results

        if energy:

            self.calc_atom_wsas()
            return self.wsas

        else:

            return self.areas

    @staticmethod
    def calc_wsas_atom_contribution(atom_series, params):
        """
        Calculate the contribution for a given atom to the estimation of the
        configurational entropy that would be computed via normal mode analysis.
        Formula from:
            Wang, J., & Hou, T. (2012). Develop and Test a Solvent Accessible
            Surface Area-Based Model in Conformational Entropy Calculations.
            Journal of Chemical Information and Modeling,
            52(5), 1199–1212. http://doi.org/10.1021/ci300064d

        contribution = weight * (sas + k * bsa)
        bsa = (4 * pi * (radius + probe_radius)^2) - sas

        sas = solvent accessible surface
        bsa = buried surface area

        Note: Units are cal/mol (NOT kcal/mol)

        Parameters
        ----------
        atom_series : pandas.core.series.Series
            Series containing 'residue', 'atom_name', 'atom_type' information
            and then data columns containing surface areas.
        params : dict
            Dictionary containing parameters for sas to energy conversion
            (originally normal mode configurational entropy estimate).

        Returns
        -------
        float
            Atomic contribution to the configurational entropy

        """

        idx = atom_series.index
        value_idx = [x for x in idx if x not in ['resid', 'residue_name',
                                                 'atom_name', 'atom_type']]
        sas = atom_series[value_idx]

        atom_type = atom_series['atom_type']

        probe_radius = params['rprobe']
        atom_type_info = params['params'][atom_type]
        weight = atom_type_info['weight']
        k = params['k']

        total_atom_surf = 4 * np.pi * (atom_type_info['radius'] + probe_radius) ** 2
        bsa = total_atom_surf - sas

        atom_series[value_idx] = weight * (sas + (k * bsa))

        return atom_series

    def calc_atom_wsas(self):


        #TODO: Check that areas exists
        areas = self.areas
        parameters = self.parameters

        for component, data in areas.items():

            self.wsas[component] = data.apply(lambda x: self.calc_wsas_atom_contribution(x, parameters), axis=1)

    @staticmethod
    def create_exclusion_filter(residues):
        """
        Convert input list of residue names into an mdtraj filter which will
        select all residues not matching an entry in the input list.

        Parameters
        ----------
        residues : list
            Residue codes to be combined into a exclusion filter.

        Returns
        -------
        str
            Selection text in mdtraj format selecting everything that does not
            match input residue names.

        """

        # TODO: Check the escape character for + ions

        selection_text = '! ('
        selection_text += ' or '.join(['(resname =~ {:s})'.format(x) for x in residues])
        selection_text += ')'

        return selection_text

    def create_component_selections(self, trajectory):
        """
        Obtain atom level selections for components used in calculation:
        complex, receptor and ligand. For component only calculations only
        one selection is returned, for complex all three.

        Parameters
        ----------
        trajectory : mdtraj.Trajectory
            An input MD trajectory - contains topology information

        Returns
        -------
        dict
            Keys are component names and values are numpy arrays of selected atom indexes
        """

        component = self.component.value

        solvent_residues = self.solvent_residues

        filters = {}

        if component in ['complex', 'ligand']:

            filters['ligand'] = self.ligand_filter

        else:

            filters['receptor'] = self.create_exclusion_filter(solvent_residues)

        if component == 'complex':
            filters['complex'] = self.create_exclusion_filter(solvent_residues)
            filters['receptor'] = '{} and (not {})'.format(filters['complex'],
                                                           filters['ligand'])

        selections = {}

        for component, filter in filters.items():
            selections[component] = trajectory.top.select(filter)

        return selections

    def update_sasa_config(self):

        parameters = self.parameters['params']
        freesasa_config_file = self.freesasa_config_file
        tmp_dir = self.tmp_dir
        nonstandard_residue_files = self.nonstandard_residue_files
        ligand_topology = self.ligand_topology

        if ligand_topology:

            files_to_add = [ligand_topology] + nonstandard_residue_files

        else:

            files_to_add = nonstandard_residue_files

        residues_to_add = {}

        for filename in files_to_add:
            residues_to_add.update(extract_residue(filename))

        if residues_to_add:

            sasa_config = os.path.join(tmp_dir, 'system_sasa.config')

            add_residues_freesasa_config_file(residues_to_add,
                                              sasa_config,
                                              parameters,
                                              orig_filename=freesasa_config_file)

            self.freesasa_config_file = sasa_config

    def compute_nm_energy(self):

        #TODO: Check have wsas

        energies = self.energies
        temperature = self.temperature
        intercept = self.parameters['intercept']

        for component, data in self.wsas.items():

            values = data.drop(['residue_name', 'resid', 'atom_name', 'atom_type'], axis =1)

            # Factor of 1000 converts cal/mol to kcal/mol
            energies[component] = temperature / 1000 * (values.sum() + intercept)

        if self.component.value == 'complex':
            energies['difference'] = energies['complex'] - energies['receptor'] - energies['ligand']


def validate_prmtop_filename(original_prmtop, target_dir='', force_new=False):
    """
    Check that file exists and create a symlink if it doesn't have a
    prmtop extension (often *.top is used but mdtraj cant't detect type
    with ambiguous extensions).

    Parameters
    ----------
    original_filename : str
        Path to supposed prmtop file
    target_dir : str
        Directory in which to create symlink if required

    Returns
    -------
    Path
        Location of verified prmtop (with potentially edited filename)

    """

    original_filename = Path(original_prmtop)

    if not os.path.isfile(original_filename):
        raise IOError()

    basename, ext = os.path.splitext(os.path.basename(original_filename))

    if ext != 'prmtop':

        top_filename = Path(os.path.join(target_dir, basename + '.prmtop'))

        if os.path.islink(top_filename) and force_new:
            os.unlink(top_filename)

        os.symlink(original_filename.absolute(), top_filename.absolute())

    else:

        top_filename = original_prmtop

    return str(top_filename)


