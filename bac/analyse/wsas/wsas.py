import shutil
import os
import argparse
import tempfile
import json

from pathlib import Path

import pandas as pd
import mdtraj
import parmed as pmd

import sasa_analysis
import freesasa_utils

import extract_residues


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


def extant_file(x):
    """
    'Type' for argparse - checks that file exists but does not open.

    Parameters
    ----------
    x : str
        Candidate file path

    Returns
    -------
    Path
        Validated path

    """

    if not os.path.isfile(x):
        # ArgumentTypeError gives a rejection message of the form:
        # error: argument input: <passed error message>
        if os.path.exists(x):
            raise argparse.ArgumentTypeError("{0} is not a file".format(x))
        else:
            raise argparse.ArgumentTypeError("{0} does not exist".format(x))

    return Path(x)


def str2bool(x):
    """
    'Type' for argparse - checks that file exists but does not open.

    Parameters
    ----------
    x : str
        Commandline choice

    Returns
    -------
    Bool
        Converted from string

    """

    if x.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif x.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


def check_prmtop(top_filename):
    """
    'Type' for argparse - checks that file exists and is a valif Amber prmtop

    Parameters
    ----------
    top_filename : str
        Candidate file path

    Returns
    -------
    Path
        Validated path

    """

    # Check file exists (and convert string to Path object)
    top_path = extant_file(top_filename)

    try:
        # Use parmed to check file is a valid Amber topology
        pmd.amber.AmberParm(top_filename)
    except:
        raise argparse.ArgumentTypeError('{:s} is not a valid prmtop file.'.format(top_filename))

    return top_path


def commandline_parser():

    #TODO: Allow selection of output directory

    defaults_dirname = os.path.dirname(os.path.realpath(__file__))
    default_config = os.path.join(defaults_dirname, 'template-input.json')
    default_wsas_filename = os.path.join(defaults_dirname, 'amber_config.txt')
    default_params_filename = os.path.join(defaults_dirname,
                                           'wsas-params-wang2012.json')

    parser = argparse.ArgumentParser(description='WSAS Calculator: Computes surface area '
                                                 'related free energy components')

    parser.add_argument('-i', dest='config_filename', default=default_config,
                        help='Input JSON configuration file', metavar='FILE',
                        type=extant_file)

    parser.add_argument('-st', dest='system_topology', required=True,
                        help='Topology for the full (usually solvated) system',
                        metavar='FILE', type=check_prmtop)

    parser.add_argument("-lt", dest='ligand_topology', required=False,
                        help='Topology for the ligand (in vaccuo).',
                        metavar='FILE', type=check_prmtop)

    parser.add_argument('-t', dest='trajectories', required=False, nargs = '+',
                        help='Trajectories of coordinates to use in analysis.',
                        metavar='FILE', type=extant_file)

    parser.add_argument('-w', dest='wsas_params', default=default_wsas_filename,
                        help='JSON file containing atom/residue parameters for WSAS analysis.',
                        metavar='FILE', type=extant_file)

    parser.add_argument('-f', dest='freesasa_params', default=default_params_filename,
                        help='Freesasa atom and residue parameters file.',
                        metavar='FILE', type=extant_file)

    parser.add_argument('-c', dest='component', choices=['complex', 'receptor', 'ligand'],
                        default='complex', help='Contents to analyse - complex/receptor/ligand.')

    parser.add_argument('-o', dest='output_dir', required=False, default='.',
                        help='Directory to store output',
                        metavar='FILE', type=Path)

    args = parser.parse_args()

    return args


def validate_prmtop_filename(original_filename, target_dir=Path('.')):
    """
    Check that file exists and create a symlink if it doesn't have a
    prmtop extension (often *.top is used but mdtraj cant't detect type
    with ambiguous extensions).

    Parameters
    ----------
    original_filename : str
        Path to supposed prmtop file
    target_dir
        Directory in which to create symlink if required

    Returns
    -------
    Path
        Location of verfified prmtop (with potentially edited filename)

    """

    if not os.path.isfile(original_filename):
        raise IOError()

    basename, ext = os.path.splitext(os.path.basename(original_filename))

    if ext != 'prmtop':

        top_filename = Path(os.path.join(target_dir, basename + '.prmtop'))
        os.symlink(original_filename.absolute(), top_filename.absolute())

    else:

        top_filename = Path(original_filename)

    return top_filename


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


def create_component_selections(traj, setup):
    """
    Obtain atom level selections for components used in calculation:
    complex, receptor and ligand. For component only calculations only
    one selection is returned, for complex all three.

    Parameters
    ----------
    traj : mdtraj.Trajectory
        An input MD trajectory - contains topology information
    setup : argparse.Namespace
        Contains user selections

    Returns
    -------
    dict
        Keys are component names and values are numpy arrays of selected atom indexes
    """

    solvent_residues = setup.solvent_residues

    filters = {}

    if setup.component in ['complex', 'ligand']:

        filters['ligand'] = setup.ligand_filter

    else:
        filters['receptor'] = parse_filter_residues(solvent_residues)

    if setup.component == 'complex':

        filters['complex'] = parse_filter_residues(solvent_residues)
        filters['receptor'] = '{} and (not {})'.format(filters['complex'],
                                                       filters['ligand'])

    selections = {}

    for component, filter in filters.items():

        selections[component] = traj.top.select(filter)

    return selections


def update_sasa_config(setup):

    if setup.ligand_topology:

        files_to_add = [setup.ligand_topology] + setup.nonstandard_residue_files

    else:

        files_to_add = setup.nonstandard_residue_files

    residues_to_add = {}

    for filename in files_to_add:
        residues_to_add.update(extract_residues.extract_residue(filename))

    if residues_to_add:

        sasa_config = setup.tmp_dir.joinpath('system_sasa.config')

        freesasa_utils.add_residues_freesasa_config_file(residues_to_add,
                                                         sasa_config)

    else:

        sasa_config = freesasa_utils.default_config_filename

    return sasa_config


def wsas_calc(setup, sasa_nm_params):

    tmp_dir = setup.tmp_dir
    trajectories = setup.trajectories
    system_topology = setup.system_topology
    stride = setup.stride
    output_dir = setup.output_dir

    output_dir.mkdir(exist_ok=True, parents=True)

    system_amber_top = pmd.load_file(str(system_topology))

    for idx, trajectory_filename in enumerate(trajectories):

        # Use stride here then hack to select first/last frames as
        # not supported simply by mdtraj (and want to minimize memory usage)
        traj = mdtraj.load(str(trajectory_filename), top=str(system_topology), stride=stride)

        n_frames = traj.n_frames
        n_orig_frames = n_frames * stride

        first_frame_input = setup.first_frame

        if setup.last_frame == -1:
            last_frame_input = n_orig_frames
        else:
            last_frame_input = setup.last_frame

        # Map frame numbers in original file to those read in (i.e. account for stride)
        frames_to_use = [x for x,y in enumerate(range(0, n_orig_frames, stride))
                         if y in range(first_frame_input, last_frame_input)]

        first_frame_strided = frames_to_use[0]
        last_frame_strided = frames_to_use[-1]

        if idx == 0:

            sasa_config = update_sasa_config(setup)

            sasa_calculator = freesasa_utils.FreesasaRunner(config=sasa_config)

            atom_selections = create_component_selections(traj, setup)

            results = {}

            for component, atom_list in atom_selections.items():
                results[component] = pd.DataFrame()
                results[component]['residue'] = [traj.top.atom(x).residue.name for x in atom_list]
                results[component]['atom_name'] = [system_amber_top.parm_data['ATOM_NAME'][x] for x in atom_list]
                results[component]['atom_type'] = [system_amber_top.parm_data['AMBER_ATOM_TYPE'][x] for x in atom_list]

        for component, atom_list in atom_selections.items():

            traj_pdb_filename = os.path.join(tmp_dir, component + '-traj.pdb')
            comp_traj = traj[first_frame_strided:last_frame_strided].atom_slice(atom_list)
            comp_traj.save(traj_pdb_filename)

            atom_areas = pd.DataFrame(sasa_calculator.run(traj_pdb_filename))

            avg_areas = atom_areas.mean()

            nm_contrib = [sasa_analysis.atom_contribution_nm(results[component]['atom_type'][ndx], sasa, sasa_nm_params) for
                          ndx, sasa in enumerate(avg_areas)]

            results[component][trajectory_filename] = nm_contrib


    return results


def read_config_file(setup):

    setup_vars = vars(setup)

    with open(setup.config_filename, 'r') as f:
        calc_config = json.load(f)

    for key, value in calc_config.items():

        if key not in setup_vars:

            setup_vars[key] = value

    return


def setup_workspace(setup):

    setup.tmp_dir = Path(tempfile.mkdtemp())

    # mdtraj recognizes filetype from extension
    # Traditionally BAC used .top instead of prmtop which breaks this
    setup.system_topology = validate_prmtop_filename(setup.system_topology,
                                                     setup.tmp_dir)

    if setup.ligand_topology:

        setup.ligand_topology = validate_prmtop_filename(setup.ligand_topology,
                                                         setup.tmp_dir)

        lig_res = extract_residues.extract_residue(setup.ligand_topology)

        setup.ligand_filter = 'resname {:s}'.format(list(lig_res.keys())[0])

    read_config_file(setup)

    #TODO: Validate contents of setup


    return


def teardown_workspace(setup):

    shutil.rmtree(setup.tmp_dir)

    return


def get_wsas_component_output(setup, sasa_nm_params, results, component):

    output_file = setup.output_dir.joinpath('{}-wsas-atom-breakdown.dat'.format(component))
    results.to_csv(output_file, sep='\t', index=False)

    wsas_component = sasa_analysis.nm_component_calc(results.select_dtypes(exclude=['object']),
                                                     setup.temperature,
                                                     sasa_nm_params['intercept'])

    output_file = setup.output_dir.joinpath('{}-wsas-summary.dat'.format(component))
    wsas_component.to_csv(output_file, sep='\t')

    return wsas_component


if __name__ == "__main__":
    # execute only if run as a script

    setup = commandline_parser()

    setup_workspace(setup)

    nm_param_filename = setup.freesasa_params

    with open(nm_param_filename, 'r') as f:
        sasa_nm_params = json.load(f)

    results = wsas_calc(setup, sasa_nm_params)

    if 'complex' in results:

        for component in ['complex', 'receptor', 'ligand']:

            wsas_component = get_wsas_component_output(setup, sasa_nm_params, results[component], component)

            if component == 'complex':

                wsas = wsas_component

            else:

                wsas += -wsas_component

        output_file = setup.output_dir.joinpath('{}-wsas-result.dat'.format(component))

        with open(output_file, 'w') as outfile:
            print('{:.3f}\t{:.3f}'.format(wsas.mean(), wsas.std()), file=outfile)

    else:

        component = list(results.keys())[0]
        wsas_component = get_wsas_component_output(setup, sasa_nm_params, results[component], component)

    output_file = setup.output_dir.joinpath('{}-wsas-summary.dat'.format(component))
    wsas_component.to_csv(output_file, sep='\t')


    teardown_workspace(setup)