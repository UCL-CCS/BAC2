import os
import parmed as pmd
from bac.builder.utils.header import HeaderInfo


class SourceStructure(object):

    def __init__(self, structure_filename):

        if os.path.splitext(structure_filename)[1] != '.pdb':
            raise NotImplementedError("We can only process PDB files at present")

        self.structure = pmd.load_file(structure_filename)
        self.header = HeaderInfo(structure_filename)
        self.chains = []
        self.chain_types = {}


