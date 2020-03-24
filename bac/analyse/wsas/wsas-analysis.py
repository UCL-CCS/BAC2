import os, sys
import shutil
import argparse
from pathlib import Path

from bac.analyse.wsas.input_checks import extant_file, check_prmtop
from bac.analyse.wsas.wsas import Wsas, validate_prmtop, DEFAULT_CONFIG_FILENAME, DEFAULT_PARAMETER_FILENAME


def commandline_parser():

    # TODO: Allow user to set solvent filter
    # TODO: Create option to use a parameters YAML file rather than commandline

    parser = argparse.ArgumentParser(description='WSAS Calculator: Computes surface area '
                                                 'related free energy components')

    parser.add_argument('-st', dest='system_topology', required=True,
                        help='Topology for the full (usually solvated) system',
                        metavar='FILE', type=check_prmtop)

    parser.add_argument('-lt', dest='ligand_topology', required=False,
                        help='Topology for the ligand (in vaccuo).',
                        metavar='FILE', type=check_prmtop)

    parser.add_argument('-lf', dest='ligand_filter', required=False,
                        help='Filter text (mdtraj DSL) to identify ligand',
                        type=str)
    parser.add_argument('-lfr', dest='ligand_filter_resname', required=False,
                        help='Resname of the ligand (mdtraj DSL) to identify ligand. e.g. "GLY"',
                        type=str)
    parser.add_argument('-t', dest='trajectories', required=False, nargs='+',
                        help='Trajectories of coordinates to use in analysis.',
                        metavar='FILE', type=extant_file)

    parser.add_argument('-te', dest='temp', default=300.0,
                        help='Temperature for energy computations',
                        type=float)

    parser.add_argument('-n', dest='non_standard_residues', default=[], nargs='+',
                        help='Topology file for non-standard residues',
                        metavar='FILE', type=extant_file)

    parser.add_argument('-w', dest='wsas_params', default=DEFAULT_PARAMETER_FILENAME,
                        help='JSON file containing atom/residue parameters for WSAS analysis.',
                        metavar='FILE', type=extant_file)

    parser.add_argument('-f', dest='freesasa_config', default=DEFAULT_CONFIG_FILENAME,
                        help='Freesasa atom and residue parameters file.',
                        metavar='FILE', type=extant_file)

    parser.add_argument('-c', dest='component', choices=['complex', 'receptor', 'ligand'],
                        default='complex', help='Contents to analyse - complex/receptor/ligand.')

    parser.add_argument('-o', dest='output_dir', required=False, default='.',
                        help='Directory to store output',
                        metavar='FILE', type=Path)

    args = parser.parse_args()

    return args


def save_outputs(wsas, output_dir):

    for component, wsas_data in wsas.wsas.items():

        output_file = os.path.join(output_dir, '{}-wsas-atom-breakdown.dat'.format(component))
        wsas_data.to_csv(output_file, sep='\t', index=False)

    for component, energy_data in wsas.energies.items():
        output_file = os.path.join(output_dir, '{}-wsas-summary.dat'.format(component))
        energy_data.to_csv(output_file, sep='\t')

    if 'complex' in wsas.energies:

        nm_approx = 0.0

        for component in ['complex', 'receptor', 'ligand']:

            energies = wsas.energies[component]

            if component == 'complex':
                nm_approx = energies
            else:
                nm_approx -= energies

        output_file = os.path.join(output_dir, 'nm-approx-summary.dat')
        nm_approx.to_csv(output_file, sep='\t')

        output_file = os.path.join(output_dir, 'nm-approx-avg.dat')
        print('{:.3f}\t{:.3f}'.format(nm_approx.mean(), nm_approx.std()), file=open(output_file, 'w'))


if __name__ == "__main__":
    # execute only if run as a script

    args = commandline_parser()

    top = args.system_topology
    lig_top = args.ligand_topology
    lig_filter = args.ligand_filter
    lig_filter_resname = args.ligand_filter_resname
    trajectories = args.trajectories
    output_dir = args.output_dir
    temperature = args.temp
    component = args.component
    parameter_file = args.wsas_params
    config_file = args.freesasa_config
    non_standard_residues = args.non_standard_residues

    top = validate_prmtop(top, output_dir, override=True)
    if lig_top is not None:
        lig_top = validate_prmtop(lig_top, output_dir, override=True)

    if component in ['complex', 'ligand']:

        if not lig_filter and not lig_top:

            print('For complex/ligand calculations a ligand filter or '
                  'topology must be supplied')
            sys.exit(1)

    wsas_calculator = Wsas(topology=top,
                           component=component,
                           ligand_topology=lig_top,
                           ligand_filter_resname=lig_filter_resname,
                           trajectories=trajectories,
                           stride=10,
                           ligand_filter=lig_filter,
                           nonstandard_residue_files=non_standard_residues,
                           temperature=temperature,
                           parameter_file=parameter_file,
                           config_file=config_file,
                           )

    wsas_calculator.calc_surface_areas()
    wsas_calculator.compute_nm_energy()

    save_outputs(wsas_calculator, output_dir)

    shutil.rmtree(wsas_calculator.tmp_dir)
