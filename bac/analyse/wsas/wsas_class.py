from pathlib import Path
from typing import List
from enum import Enum
import tempfile

import parmed as pmd
import mdtraj as md
import pandas as pd

from bac.utils.decorators import advanced_property, pathlike
from .extract_residues import extract_residue
from .sasa_analysis import atom_contribution_nm
from .freesasa_utils import FreesasaRunner, add_residues_freesasa_config_file, default_config_filename


class Component(Enum):
    complex = 'complex'
    ligand = 'ligand'
    receptor = 'receptor'


class Wsas:
    def __init__(self, **kwargs):
        self.temperature = kwargs.get('temperature')
        self.slice = slice(kwargs.get('first_frame', 0), kwargs.get('last_frame', -1), kwargs.get('stride', 1))
        self.ligand_filter = kwargs.get('ligand_filter')
        self.nonstandard_residue_files = kwargs.get('nonstandard_residue_files')
        self.solvent_residues = kwargs.get('solvent_residues')

        self.topology = kwargs.get('topology')
        self.ligand_topology = kwargs.get('ligand_topology')
        self.trajectories: List[md.Trajectory] = kwargs.get('trajectories')
        self.component = kwargs.get('component')
        self.output = kwargs.get('output')
        self.tmp_dir = tempfile.mkdtemp()

    @advanced_property(type=pmd.amber.AmberParm)
    def topology(self): pass

    @advanced_property(type=pmd.amber.AmberParm)
    def ligand_topology(self): pass

    @ligand_topology.post_set_processing
    def ligand_topology(self):
        # Turn this file into a list of residue name: atom names.
        # Then set the ligand filter to be this. lol.
        pass

    @advanced_property(type=Component, default=Component.complex)
    def component(self): pass

    @advanced_property(type=str)
    def ligand_filter(self):
        # Read default from the configuration file pass somehow.
        pass

    @pathlike(default='.')
    def output(self): pass

    def calculate_wsas(self, sasa_nm_params):
        self.output.mkdir(exist_ok=True, parents=True)

        any_trajectory = self.trajectories[0]

        sasa_config = self.update_sasa_config()
        sasa_calculator = FreesasaRunner(config=sasa_config)
        atom_selections = self.create_component_selections(any_trajectory)

        results = {}

        for component, atom_list in atom_selections.items():
            results[component] = pd.DataFrame()
            results[component]['residue'] = [any_trajectory.top.atom(x).residue.name for x in atom_list]
            results[component]['atom_name'] = [self.topology.parm_data['ATOM_NAME'][x] for x in atom_list]
            results[component]['atom_type'] = [self.topology.parm_data['AMBER_ATOM_TYPE'][x] for x in atom_list]

        for trajectory in self.trajectories:

            for component, atom_list in atom_selections.items():

                trajectory_pdb = '{}/{}-traj.pdb'.format(self.tmp_dir, component)
                trajectory[self.slice].save_pdb(trajectory_pdb)

                atom_areas = pd.DataFrame(sasa_calculator.run(trajectory_pdb))

                avg_areas = atom_areas.mean()

                nm_contrib = [
                    atom_contribution_nm(results[component]['atom_type'][ndx], sasa, sasa_nm_params) for
                    ndx, sasa in enumerate(avg_areas)]

                results[component][repr(trajectory)] = nm_contrib

        return results

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

        filters = {}

        if self.component in ['complex', 'ligand']:
            filters['ligand'] = self.ligand_filter

        else:
            filters['receptor'] = self.parse_filter_residues(self.solvent_residues)

        if self.component == 'complex':
            filters['complex'] = self.parse_filter_residues(self.solvent_residues)
            filters['receptor'] = '{} and (not {})'.format(filters['complex'],
                                                           filters['ligand'])

        selections = {}

        for component, filter in filters.items():
            selections[component] = trajectory.top.select(filter)

        return selections

    def update_sasa_config(self):

        if self.ligand_topology:

            files_to_add = [self.ligand_topology] + self.nonstandard_residue_files

        else:

            files_to_add = self.nonstandard_residue_files

        residues_to_add = {}

        for filename in files_to_add:
            residues_to_add.update(extract_residue(filename))

        if residues_to_add:

            sasa_config = self.tmp_dir.joinpath('system_sasa.config')

            add_residues_freesasa_config_file(residues_to_add, sasa_config)

        else:

            sasa_config = default_config_filename

        return sasa_config

    @staticmethod
    def parse_filter_residues(filter_residues):
        """
        Convert input list of residue names into an mdtraj filter which will
        select all residues not matching an entry in the input list.

        Parameters
        ----------
        filter_residues : list
            Residue codes to be combined into a exclusion filter.

        Returns
        -------
        str
            Selection text in mdtraj format selecting everything that does not
            match input residue names.

        """

        # TODO: Check the escape character for + ions

        selection_text = '! ('
        selection_text += ' or '.join(['(resname =~ {:s})'.format(x) for x in filter_residues])
        selection_text += ')'

        return selection_text



