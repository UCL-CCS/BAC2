import os
import json
import tempfile
from enum import Enum
from pkg_resources import resource_filename

import numpy as np
import mdtraj as md
import pandas as pd
import parmed as pmd

from bac.utils.decorators import advanced_property
from .extract_residues import extract_residue
from .freesasa_utils import FreesasaRunner


class Component(Enum):
    complex = 'complex'
    ligand = 'ligand'
    receptor = 'receptor'

# Defaults


DEFAULT_PARAMETER_FILENAME = resource_filename(__name__, 'data/wsas-params-wang2012.json')
DEFAULT_CONFIG_FILENAME = resource_filename(__name__, 'data/amber_config.txt')


class Wsas:

    def __init__(self, component, trajectories, topology, ligand_topology=None, ligand_filter=None, temperature=300,
                 first_frame=0, last_frame=None, stride=1, nonstandard_residue_files=[], solvent_residues=None,
                 parameter_file=None, config_file=None):
        """Wsas class used to analyse simulations.

        Parameters
        ----------
        topology: str
            File path to the complex topology file.
        trajectories: list of str
            List of paths to trajectory files.
        component: Component
            Which component to analyse. One of complex, ligand or receptor.
        temperature: float
            Temperature of the system in kelvin. Default is 300K.
        first_frame: int
            Index of the first frame of every trajectory to analyse. Default is the 0, i.e. first frame.
        last_frame: int
            Index of the last frame of every trajectory to analyse. Default is -1, i.e. the last frame.
        stride: int
            Stride of frame reading. Use this if you want to skip every `stride` frames.
        ligand_filter
        nonstandard_residue_files
        solvent_residues
        ligand_topology
        parameter_file
        config_file
        """

        self.temperature = temperature
        self.slice = slice(first_frame, last_frame, stride)
        self.ligand_filter = ligand_filter
        self.nonstandard_residue_files = nonstandard_residue_files

        default_solvent = ["WAT", "HOH", "'Cl.*'", "CIO", "'Cs+'", "IB", "'K.*'",
                           "'Li+'", "'MG.*'", "'Na+'", "'Rb+'", "CS", "RB", "NA",
                           "F",  "ZN"]
        self.solvent_residues = solvent_residues or default_solvent

        self.topology = topology
        self.ligand_topology = ligand_topology
        self.trajectories = trajectories
        self.component = component

        self.tmp_dir = tempfile.mkdtemp()

        self.areas = {}
        self.wsas = {}
        self.energies = {}

        with open(parameter_file or DEFAULT_PARAMETER_FILENAME) as f:
            self.parameters = json.load(f)

        self.freesasa_config_file = config_file or DEFAULT_CONFIG_FILENAME

        # self.ligand_filter = self.ligand_filter or 'resname {:s}'.format(
        # list(extract_residue(self.ligand_topology).keys())[0]) if self.ligand_topology else None

        if not self.ligand_filter and self.ligand_topology:

            lig_res = extract_residue(self.ligand_topology)

            self.ligand_filter = 'resname {:s}'.format(list(lig_res.keys())[0])

    @advanced_property(type=str)
    def topology(self): pass

    @advanced_property(type=Component, default=Component.complex)
    def component(self): pass

    def calc_surface_areas(self, wsas=True):
        """
        Calculate solvent accessible surface area (SASA), and optionally derived
        weighted surface areas, for trajectories in self.trajectories.

        Notes
        -----
        Outputs are placed in a dict of pd.DataFrames where the key is the
        component ('complex'/'receptor'/'ligand'). SASA results are stored in
        self.areas, weighted in self.wsas.

        Parameters
        ----------
        wsas : bool
            Should energy estimates be computed?

        Returns
        -------
        dict
            Final results - either self.areas or self.wsas depending on
            whether the later are calculated

        """

        topology_filename = self.topology
        top = pmd.load_file(topology_filename)

        sasa_calculator = FreesasaRunner(config=self.freesasa_config_file, ligand_topology=self.ligand_topology,
                                         wsas_params=self.parameters['params'], tmp_dir=self.tmp_dir,
                                         nonstandard_residue_files=self.nonstandard_residue_files)

        atom_selections = None

        results = {}

        for trajectory in self.trajectories:

            traj = md.load(trajectory, top=topology_filename)

            # If this is the first trajectory setup data frames to store data
            if atom_selections is None:

                atom_selections = self.create_component_selections(traj)

                for component, atom_list in atom_selections.items():
                    results[component] = pd.DataFrame()
                    results[component]['resid'] = [top.atoms[x].residue.number for x in atom_list]
                    results[component]['residue_name'] = [top.atoms[x].residue.name for x in atom_list]
                    results[component]['atom_name'] = [top.atoms[x].name for x in atom_list]
                    missing_params = []
                    for x in atom_list:
                        if top.atoms[x].type not in self.parameters['params']:
                            missing_params.append(top.atoms[x])

                    if missing_params:
                        raise Exception('wsas calculation will not work because there is no such type', ', '.join(missing_params))
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

        if wsas:

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
        """
        Calculate weighted SASA (wsas) values from surface area for each
        component with SASA values in self.areas using
        self.calc_wsas_atom_contribution (see docstring there for details of
        calculation).

        Parameters
        ----------

        Notes
        -----
        Updates self.wsas

        Returns
        -------

        """

        # TODO: Check that areas exists
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

    def compute_nm_energy(self):
        """
        Compute estimates of the energy associated with normal mode
        configurational entropy from atomic wsas values for each component.
        For 'complex' compute difference; complex - receptor -ligand
        Calculation is:

        energy = T * S_nm = T * (sum(wsas) + intercept)

        Intercept is read in from the input wsas parameters (via
        self.parameters)

        Parameters
        ----------

        Notes
        -----
        Updates self.energies with pd.Dataframes for each component (plus
        difference for complex component)

        Returns
        -------

        """

        # TODO: Check have wsas

        energies = self.energies
        temperature = self.temperature
        intercept = self.parameters['intercept']

        for component, data in self.wsas.items():

            values = data.drop(['residue_name', 'resid', 'atom_name', 'atom_type'], axis=1)

            # Factor of 1000 converts cal/mol to kcal/mol
            energies[component] = temperature / 1000 * (values.sum() + intercept)

        if self.component.value == 'complex':
            energies['difference'] = energies['complex'] - energies['receptor'] - energies['ligand']


def validate_prmtop(prmtop, target_dir=None, override=False):
    """
    Check that file exists and create a symlink if it doesn't have a
    prmtop extension (often *.top is used but mdtraj cant't detect type
    with ambiguous extensions).

    Parameters
    ----------
    prmtop : str
        Path to supposed prmtop file
    target_dir : str
        Directory in which to create symlink if required. If None the symlink will
        be created in the same directory.
    override: bool
        Override possible already existing file with .prmtop extension. Default is false

    Returns
    -------
    Path
        Location of verified prmtop (with potentially edited filename)
    """

    if not os.path.isfile(prmtop):
        raise IOError()

    _, ext = os.path.splitext(prmtop)

    if ext is 'prmtop':
        return prmtop

    target_dir = target_dir or os.path.dirname(os.path.abspath(prmtop))

    new_prmtop = os.path.join(target_dir, os.path.basename(prmtop) + '.prmtop')

    if os.path.islink(new_prmtop) and override:
        os.unlink(new_prmtop)

    os.symlink(os.path.abspath(prmtop), os.path.abspath(new_prmtop))

    return new_prmtop
