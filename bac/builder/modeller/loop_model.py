basic_modeller_script = """
from modeller import *
from modeller.automodel import *    # Load the automodel class

log.verbose()
env = environ()

class MyModel(automodel):
    def select_atoms(self):
        return selection(self.residue_range('133', '135'),
                         self.residue_range('217', '230'))

a = MyModel(env, alnfile = 'alignment.ali',
            knowns = '1qg8', sequence = '1qg8_fill')
a.starting_model= 1
a.ending_model  = 1

a.make()
"""