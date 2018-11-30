import os
from collections import OrderedDict
from parmed.amber.offlib import AmberOFFLibrary as off_parser

# TODO: Make this sane
amberhome = os.environ['AMBERHOME']
leaprc_path = os.path.join(amberhome, 'dat', 'leap', 'cmd')


def leaprc_off_files(leaprcs):

    off_files = []

    for leaprc in leaprcs:
        leaprc_file = os.path.join(leaprc_path, leaprc)
        with open(leaprc_file) as rc:
            for line in rc:
                if line.startswith('loadOff'):
                    off_files.append(line.split()[1])

    return off_files


def off_to_residues(off_files):

    ff_residues = OrderedDict()

    for off_file in off_files:

        ff_residues.update(off_parser.parse(off_file))

    return ff_residues
