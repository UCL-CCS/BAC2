import os
import parmed as pmd
from bac.builder.utils.header import HeaderInfo
from bac.builder.structure_utils import scan_chain_type, update_chain_type_assignment


class SourceStructure(object):

    def __init__(self, structure_filename):

        if os.path.splitext(structure_filename)[1] != '.pdb':
            raise NotImplementedError("We can only process PDB files at present")

        self.structure = pmd.load_file(structure_filename)
        self.header = HeaderInfo(structure_filename)
        self.chains = set([x.chain for x in self.structure.residues])
        self.chain_types = {}

        chain_types = self.chain_types

        for chain in self.chains:
            chain_types[chain] = scan_chain_type(self, chain)
            update_chain_type_assignment(self, chain_types[chain])


