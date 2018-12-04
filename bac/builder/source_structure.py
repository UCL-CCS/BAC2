import os
import numpy as np
import parmed as pmd
from bac.builder.utils.header import HeaderInfo
from bac.builder.structure_utils import (scan_chain_type, update_chain_type_assignment,
                                         get_polymer_gaps)
from bac.builder.utils.sequence import convert_resname_list
import re


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
        structure = self.structure

        for chain in self.chains:
            chain_types[chain] = scan_chain_type(structure, chain)

            # Record what is nonstandard before we update types to consider
            # bonding for residues linked to polymers (protein, RNA, DNA).
            condition = (chain_types[chain][:, 1] == 'unknown')
            self.chain_nonstandard[chain] = chain_types[chain][condition][:, 0]

            # Update types to classify nonstandard residues bonding in
            # polymers as part of neighbouring chain type.
            update_chain_type_assignment(structure, chain_types[chain])

            # Gaps defined by pair of polymer residues where residue number
            # different by > 1 and no polymer bond between them.
            self.chain_gaps[chain] = get_polymer_gaps(structure, chain,
                                                      chain_types[chain])

        self.decomposition_mapping = {}
        self.model_selection = {}

    def generate_decomposition(self):
        """
        Decompose initial chains into modelling/simulation friendly units.

        Returns
        -------

        """

        chain_types = self.chain_types
        chain_gaps = self.chain_gaps

        decomposition = {}

        for chain_id, residue_types in chain_types.items():

            decomposition[chain_id] = []

            type_changes = np.where(residue_types[:-1, 1] !=
                                    residue_types[1:, 1])[0] + 1

            type_blocks = np.split(residue_types, type_changes)

            for type_block in type_blocks:

                # Use indices to find blocks of residues with
                # contiguous numbers (indicative of a real biological chain)
                number_residue_blocks = self.residue_no_blocks(type_block[:, 0])

                for number_residue_block in number_residue_blocks:

                    types_data = residue_types[np.where(
                        np.isin(residue_types[:, 0], number_residue_block))]

                    decomposition[chain_id].append(types_data)

        return decomposition

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

    def check_residue_no_increases(self, idxs):
        """
        Check that residue with given indicies have monotonically increasing
        residue numbering.

        Parameters
        ----------
        idxs : list or np.array
            Indices of the residues to check

        Returns
        -------
        bool
            Do the relevant residue numbers increase monotonically?
        """

        residues = self.structure.residues

        residue_numbers = [residue.number for residue
                           in residues if residue.idx in idxs]

        dx = np.diff(residue_numbers)
        return np.all(dx > 0)

    def residue_no_blocks(self, idxs):
        """
        Get list of blocks of residue indexes where the residue numbers
        increase.

        Parameters
        ----------
        idxs : list or np.array
            Indices of the residues to group

        Returns
        -------
        list
            List of `np.arrays` containing indices of residues which increase in
            residue numbering.
        """

        residues = self.structure.residues

        residue_numbers = [(residue.idx, residue.number) for residue
                           in residues if residue.idx in idxs]

        residue_numbers = np.array(residue_numbers)

        return [arr[:, 0] for arr in
                np.split(residue_numbers,
                         np.where(np.diff(residue_numbers[:, 1]) < 1)[0] + 1)
                if arr.size > 1]

    def residue_contiguous_no_blocks(self, idxs):
        """
        Get list of blocks of residue indexes where the residue numbers
        are contiguous (currently assume same residue numbers are insertions).

        Parameters
        ----------
        idxs : list or np.array
            Indices of the residues to group

        Returns
        -------
        list
            List of `np.arrays` containing indices of residues which increase in
            residue numbering.
        """

        residues = self.structure.residues

        residue_numbers = [(residue.idx, residue.number) for residue
                           in residues if residue.idx in idxs]

        residue_numbers = np.array(residue_numbers)

        diff = np.diff(residue_numbers[:, 1])

        return [arr[:, 0] for arr in
                np.split(residue_numbers,
                         np.where(np.logical_or((diff < 0),
                                                (diff > 1)))[0] + 1)
                if arr.size > 1]

    def chain_sequence(self, chain_id, seq_format='letter', gap_char='-'):
        """
        Get sequence information for selected chain in input structure.

        Parameters
        ----------
        chain_id : str
            Selected chain identifier.
        seq_format : str
            Format for output - 'resname' (list of three letter reside codes),
            'letter' (list of single letter residue codes), 'fasta' (sting of
            single letter codes).
        gap_char : str
            Character to use in single letter format for gaps.

        Returns
        -------
        list or str
            Output of all residue codes in the format selected in seq_format.
        """

        if len(gap_char) > 1:
            raise RuntimeError("Must use a single character to represent "
                               "sequence gaps")

        struct = self.structure

        residues = self.structure.residues
        residue_types = self.chain_types[chain_id]

        if chain_id in self.chain_gaps:
            gaps = self.chain_gaps[chain_id]
        else:
            gaps = []

        if not np.any(np.isin(residue_types[:, 1], ['protein', 'dna', 'rna'])):
            if seq_format == 'fasta':
                return ''
            else:
                return []

        mask = np.isin(residue_types[:, 1], ['protein', 'dna', 'rna'])
        seq_idxs = residue_types[:, 0][mask]

        residue_blocks = self.residue_contiguous_no_blocks(seq_idxs)

        gap_fill = {}

        for gap in gaps:

            start = residues[gap[0]].number
            end = residues[gap[1]].number

            gap_fill[start] = [gap_char for _ in range(end - start - 1)]

        if seq_format == 'fasta':
            sequence = ''
        else:
            sequence = []

        for residue_block in residue_blocks:

            block_sequence = [residues[x].name for x in residue_block]

            if seq_format != 'resname':
                block_sequence = convert_resname_list(block_sequence,
                                                      seq_format=seq_format)

            last_residue = residues[residue_block[-1]].number

            sequence += block_sequence

            if last_residue in gap_fill:
                if seq_format == 'fasta':
                    sequence += ''.join(gap_fill[last_residue])
                else:
                    sequence += gap_fill[last_residue]

        return sequence

    def reconcile_structure_sequence(self, chain_id, sequence):

        alignment = ()

        seq_regex = self.chain_sequence(chain_id,
                                        seq_format='fasta', gap_char='.')

        matches = re.search(seq_regex, sequence)

        if matches:
            pre_seq, post_seq = re.split(seq_regex, sequence)
            alignment = (matches[0], pre_seq, post_seq)

        return alignment
