import os
import parmed as pmd
from bac.builder.utils.header import HeaderInfo
from bac.builder.structure_utils import (scan_chain_type, update_chain_type_assignment,
                                         get_chain_number_gaps)


class SourceStructure(object):

    def __init__(self, structure_filename):

        if os.path.splitext(structure_filename)[1] != '.pdb':
            raise NotImplementedError("We can only process PDB files for now.")

        self.structure = pmd.load_file(structure_filename)
        self.header = HeaderInfo(structure_filename)
        self.chains = set([x.chain for x in self.structure.residues])

        self.chain_types = {}
        self.chain_nonstandard = {}
        self.chain_gaps = {}

        chain_types = self.chain_types

        for chain in self.chains:
            chain_types[chain] = scan_chain_type(self, chain)

            # Record what is nonstandard before we update types to consider
            # bonding for residues linked to polymers (protein, RNA, DNA).
            condition = (chain_types[chain][:, 1] == 'unknown')
            self.chain_nonstandard[chain] = chain_types[chain][condition][:, 0]

            # Update types to classify nonstandard residues bonding in
            # polymers as part of neighbouring chain type.
            update_chain_type_assignment(self, chain_types[chain])

            self.chain_gaps[chain] = get_chain_number_gaps(self, chain)
