from pathlib import Path

import parmed as pmd
import mdtraj as md
import pandas as pd

from bac.utils.decorators import advanced_property


class Wsas:
    def __init__(self, **kwargs):
        self.temperature = kwargs.get('temperature')
        self.first_frame = kwargs.get('first_frame')
        self.last_frame = kwargs.get('last_frame')
        self.stride = kwargs.get('stride')
        self.ligand_filter = kwargs.get('ligand_filter')
        self.nonstandard_residue_files = kwargs.get('nonstandard_residue_files')
        self.solvent_residues = kwargs.get('solvent_residues')

        self.topology = kwargs.get('topology')
        self.ligand_topology = kwargs.get('ligand_topology')
        self.trajectories = kwargs.get('trajectories')
        self.components = kwargs.get('components')
        self.output = kwargs.get('output')

    @advanced_property(type=pmd.amber.AmberParm)
    def topology(self): pass

    @advanced_property(type=md.Topology)
    def ligand_topology(self): pass

    @ligand_topology.post_set_processing
    def ligand_topology(self):
        # Turn this file into a list of residue name: atom names.
        # Then set the ligand filter to be this. lol.
        pass

    @advanced_property(type=str)
    def ligand_filter(self):
        # Read default from the configuration file pass somehow.
        pass

    def calculate_wsas(self):
        self.output.mkdir(exist_ok=True, parents=True)

        sasa_config = update_sasa_config(setup)
        sasa_calculator = freesasa_utils.FreesasaRunner(config=sasa_config)

        results = {}

        for idx, traj in enumerate(self.trajectories):

            # Use stride here then hack to select first/last frames as
            # not supported simply by mdtraj (and want to minimize memory usage)
            # traj = mdtraj.load(str(trajectory_filename), top=str(system_topology), stride=stride)

            n_frames = traj.n_frames
            n_orig_frames = n_frames * self.stride

            first_frame_input = self.first_frame

            if self.last_frame == -1:
                last_frame_input = n_orig_frames
            else:
                last_frame_input = self.last_frame

            # Map frame numbers in original file to those read in (i.e. account for stride)
            frames_to_use = [x for x, y in enumerate(range(0, n_orig_frames, self.stride))
                             if y in range(first_frame_input, last_frame_input)]

            first_frame_strided = frames_to_use[0]
            last_frame_strided = frames_to_use[-1]

            if idx == 0:

                atom_selections = create_component_selections(traj, setup)

                for component, atom_list in atom_selections.items():
                    results[component] = pd.DataFrame()
                    results[component]['residue'] = [traj.top.atom(x).residue.name for x in atom_list]
                    results[component]['atom_name'] = [system_amber_top.parm_data['ATOM_NAME'][x] for x in atom_list]
                    results[component]['atom_type'] = [system_amber_top.parm_data['AMBER_ATOM_TYPE'][x] for x in
                                                       atom_list]

            for component, atom_list in atom_selections.items():
                traj_pdb_filename = os.path.join(tmp_dir, component + '-traj.pdb')
                comp_traj = traj[first_frame_strided:last_frame_strided].atom_slice(atom_list)
                comp_traj.save(traj_pdb_filename)

                atom_areas = pd.DataFrame(sasa_calculator.run(traj_pdb_filename))

                avg_areas = atom_areas.mean()

                nm_contrib = [
                    sasa_analysis.atom_contribution_nm(results[component]['atom_type'][ndx], sasa, sasa_nm_params) for
                    ndx, sasa in enumerate(avg_areas)]

                results[component][trajectory_filename] = nm_contrib

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

        if self.components in ['complex', 'ligand']:

            filters['ligand'] = setup.ligand_filter

        else:
            filters['receptor'] = parse_filter_residues(self.solvent_residues)

        if setup.component == 'complex':
            filters['complex'] = parse_filter_residues(self.solvent_residues)
            filters['receptor'] = '{} and (not {})'.format(filters['complex'],
                                                           filters['ligand'])

        selections = {}

        for component, filter in filters.items():
            selections[component] = trajectory.top.select(filter)

        return selections

