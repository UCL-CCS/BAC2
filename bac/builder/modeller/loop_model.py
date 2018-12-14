import os
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

a = MyModel(env, alnfile = $alignment_file,
            knowns = $known_structure, sequence = $sequence_id')
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

a = MyModel(env, alnfile = $alignment_file,
            knowns = $known_structure, sequence = $sequence_id')
a.starting_model= 1
a.ending_model  = 1

a.make()

""")


def create_loop_model_script(structure_filename, sequence_id,
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
        lines.append(" " * 25 + f"self.residue_range('{res_range[0]}', "
                                f"'{res_range[1]}')")

    sel = ',\n'.join(lines)

    script_text = template_text.substitute(known_structure=structure_name,
                                           sequence_id=sequence_id,
                                           alignment_file=alignment_file,
                                           selection=sel)

    script_name = f'{structure_name}-loop-model.py'
    script_path = os.path.join(workdir, script_name)

    with open(script_path, 'w') as fout:
        fout.write(script_text)

    return script_path
