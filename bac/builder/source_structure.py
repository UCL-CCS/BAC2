import os
import parmed as pmd
from bac.builder.utils.header import HeaderInfo
from bac.builder.structure_utils import (scan_chain_type, update_chain_type_assignment,
                                         get_polymer_gaps)


class SourceStructure(object):

    def __init__(self, structure_filename):

        if os.path.splitext(structure_filename)[1] != '.pdb':
            raise NotImplementedError("We can only process PDB files for now.")

        self.source_file = os.path.realpath(os.path.expanduser(structure_filename))
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

            # Gaps defined by pair of polymer residues where residue number
            # different by > 1 and no polymer bond between them.
            self.chain_gaps[chain] = get_polymer_gaps(self, chain,
                                                      chain_types[chain])

        self.decomposition_mapping = {}
        self.model_selection = {}

    def generate_decomposition(self):
        """
        Decompose initial chains into modelling/simulation friendly units.

        Returns
        -------

        """

        pass

    def write_decomposition(self, filename):
        """
        Write decomposition information to YAML file

        Returns
        -------

        """

        pass

    def create_scaffold(self, decomposition, selection, naming_scheme):

        pass

    def reconcile_header_sequence(self):

        pass

