import os
import textwrap
from string import Template

SELECTION_INDENT_LENGTH = 25

BASIC_MODEL_SCRIPT = Template("""
from modeller import *
from modeller.automodel import *    # Load the automodel class

log.verbose()
env = environ()

class MyModel(automodel):
    def select_atoms(self):
        # Note these are numbers of residues that change
        # Numbering starts at 1 (i.e. can be offset from PDB numbering)
        return selection(
$selection
                         )

a = MyModel(env, alnfile = '$alignment_file',
            knowns = '$known_structure', sequence = '$sequence_id')
a.starting_model= 1
a.ending_model  = 1

a.make()

""")

STANDARD_MODEL_SCRIPT = Template("""
from modeller import *
from modeller.automodel import *    # Load the automodel class

log.verbose()
env = environ()

class MyModel(loopmodel):
    def select_atoms(self):
        # Note these are numbers of residues that change
        # Numbering starts at 1 (i.e. can be offset from PDB numbering)
        return selection(
$selection
                         )

a = MyModel(env, alnfile = '$alignment_file',
            knowns = '$known_structure', sequence = '$sequence_id')
a.starting_model= 1
a.ending_model  = 1

a.make()

""")


def create_loop_model_script(structure_filename,
                             loop_residue_ranges,
                             alignment_file='alignment.ali',
                             workdir='.', basic=False):

    if structure_filename[-4:] == '.pdb':
        structure_name = structure_filename[:-4]
    else:
        structure_name = structure_filename

    if basic:
        template_text = BASIC_MODEL_SCRIPT
    else:
        # Uses modeller `loopmodel` which refines the loops
        # Needed for long loops and terminal extensions
        template_text = STANDARD_MODEL_SCRIPT

    lines = []
    for res_range in loop_residue_ranges:
        lines.append(" " * SELECTION_INDENT_LENGTH +
                     f"self.residue_range('{res_range[0]}', "
                     f"'{res_range[1]}')")

    sel = ',\n'.join(lines)

    script_text = template_text.substitute(known_structure=structure_name,
                                           sequence_id=structure_name + '_fill',
                                           alignment_file=alignment_file,
                                           selection=sel)

    script_name = f'{structure_name}-loop-model.py'
    script_path = os.path.join(workdir, script_name)

    with open(script_path, 'w') as fout:
        fout.write(script_text)

    return script_name


def write_ali_file(filename, target_seq, structure_seq,
                   structure_id='pdb', start_resid=1, end_resid=-1, chain='?'):
    with open(filename, 'w') as output_file:

        output_file.write(f'>P1;{structure_id}\n')
        output_file.write(f'structureX:{structure_id}:{start_resid}:{chain}:{end_resid}:{chain}::::\n')
        for line in textwrap.wrap(structure_seq + '*', width=75, break_on_hyphens=False):
            output_file.write(line + '\n')
        output_file.write(f'>P1;{structure_id}_fill\n')
        output_file.write('sequence:::::::::\n')
        for line in textwrap.wrap(target_seq + '*', width=75, break_on_hyphens=False):
            output_file.write(line + '\n')
