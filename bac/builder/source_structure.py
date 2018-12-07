import os
from collections import defaultdict
import pprint
import re
import parmed as pmd
from bac.builder.utils.header import HeaderInfo
from bac.builder.utils.sequence import convert_resname_list
from bac.builder.structure_utils import (get_residue_type,
                                         check_residue_polymer_bonded,
                                         align_sequence_structure,
                                         residue_list_sequence)


class Subdivision(object):

    def __init__(self, residue_type=None, chain_id='', structure=None):

        self.residues = pmd.TrackedList()
        self.residue_type = residue_type
        self.chain_id = chain_id
        self.gaps = []
        self.aligned_sequence = ''
        self.pre_sequence = ''
        self.post_sequence = ''
        self.src_structure = structure

    def is_polymer(self):

        if self.residue_type in ['protein', 'dna', 'rna']:
            return True

        return False

    def sequence(self, gap_char='-', seq_format='fasta'):

        return residue_list_sequence(self.residues, seq_format=seq_format,
                                     gap_char=gap_char)

    def align_sequence(self, sequence):

        if not self.is_polymer():
            return False

        aligned, pre, post = align_sequence_structure(self.residues, sequence)

        if not aligned:
            return False
        else:
            self.aligned_sequence = aligned
            self.pre_sequence = pre
            self.post_sequence = post
            return True

    def get_added_residue_positions(self, pre_res=0, post_res=0,
                                    numbering='modeller'):

        if numbering == 'modeller':
            offset = 0
        else:
            offset = self.residues[0].number - 1

        sequence = '-'*pre_res + self.sequence() + '-'*post_res

        return [(m.start(0) + 1 + offset, m.end(0) + offset) for m
                in re.finditer('-+', sequence)]

    def target_sequence(self, pre_res=0, post_res=0):

        aligned_sequence = self.aligned_sequence

        if aligned_sequence:

            return (self.pre_sequence[:pre_res] + aligned_sequence +
                    self.post_sequence[:post_res])

        else:
            return ''

    def modeller_selection(self):
        pass

    def __repr__(self):

        residues = self.residues
        residue_type = self.residue_type
        if self.aligned_sequence:
            aligned_sequence = self.aligned_sequence
        else:
            aligned_sequence = None

        if residues:

            first_idx = residues[0].idx
            last_idx = residues[-1].idx

        else:

            first_idx = None
            last_idx = None

        return (f"{self.__class__}( residue_idxs={first_idx}-{last_idx}," 
                f" residue_type={residue_type}, aligned_sequence={aligned_sequence})")

    def __str__(self):

        return pprint.pformat(self.__dict__, indent=4)


class SourceStructure(object):

    def __init__(self, structure, header_filename=None):

        if isinstance(structure, pmd.Structure):
            self.structure = structure
            self.source_file = ''

        else:

            if os.path.splitext(structure)[1] != '.pdb':
                raise NotImplementedError("We can only process PDB files "
                                          "for now.")

            if header_filename is None:
                header_filename = structure

            self.source_file = os.path.realpath(os.path.expanduser(structure))
            self.structure = pmd.load_file(structure)

        self.header_filename = header_filename
        self.header = HeaderInfo(header_filename)

        struct = self.structure

        self.chain_map = defaultdict(pmd.TrackedList)
        self.nonstandard = []

        for r in struct.residues:

            residue_type = get_residue_type(r)

            if residue_type == 'unknown':
                self.nonstandard.append(r)
                polymer_bond = check_residue_polymer_bonded(r)
                if polymer_bond is not None:
                    residue_type = polymer_bond

            # TODO: Is there a less hacky way to do this?
            r.residue_type = residue_type
            self.chain_map[r.chain].append(r)

        self.chains = self.chain_map.keys()
        self.decomposition = self.default_decomposition()

    def __repr__(self):

        return (f"{self.__class__}(structure={self.structure}, "
                f"source_file={self.source_file}, "
                f"header_file={self.header_filename})")

    def default_decomposition(self):

        structure = self.structure
        chain_map = self.chain_map
        decomposition = defaultdict(None)

        for chain, chain_residues in chain_map.items():

            decomposition[chain] = []
            subdivision = None

            for residue in chain_residues:

                residue_type = residue.residue_type

                if subdivision is None:

                    subdivision = Subdivision(residue_type=residue_type,
                                              chain_id=chain,
                                              structure=structure)

                elif (residue_type != subdivision.residue_type or
                      residue.number < subdivision.residues[-1].number):

                    decomposition[chain].append(subdivision)
                    subdivision = Subdivision(residue_type=residue_type,
                                              chain_id=chain,
                                              structure=structure)

                subdivision.residues.append(residue)

                if len(subdivision.residues) > 1:
                    last_residue = subdivision.residues[-2]
                    if residue.number - last_residue.number > 1:

                        bonded = check_residue_polymer_bonded(residue,
                                                              bond_type=residue_type,
                                                              directions=[-1],
                                                              specific_partner=last_residue)

                        if not bonded:
                            subdivision.gaps.append((last_residue, residue))

            if residue == chain_residues[-1]:
                decomposition[chain].append(subdivision)

        return decomposition

    def write_decomposition(self, filename):
        """
        Write decomposition information to YAML file

        Returns
        -------

        """

        pass

    def create_scaffold(self, selection, naming_scheme):

        pass

    def chain_sequence(self, chain_id, seq_format='fasta', src='structure'):

        if src not in ['structure', 'header']:
            raise RuntimeError("src must be 'structure', 'header'")

        decomposition = self.decomposition
        header_info = self.header

        if src == 'header':

            three_letter_seq = header_info.sequence

            if seq_format == 'fasta':

                sequence = convert_resname_list(three_letter_seq, seq_format=seq_format)

            elif seq_format == 'letter':

                sequence = convert_resname_list(three_letter_seq)

            else:

                sequence = three_letter_seq

            return sequence

        else:

            if seq_format == 'fasta':
                sequence = ''
            else:
                sequence = []

            if decomposition[chain_id] is None:
                return sequence

            # TODO: Warning for multi-component systems?
            if len(decomposition[chain_id]) > 1:
                end_char = '*'
            else:
                end_char = ''

            for subdivision in decomposition[chain_id]:

                sequence += residue_list_sequence(subdivision.residues,
                                                  seq_format=seq_format)
                if seq_format == 'fasta':
                    sequence += end_char

            return sequence
