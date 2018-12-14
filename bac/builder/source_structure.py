"""
The aim of the classes here is to provide a layer on top of ParmEd to help read
in structure files (currently on PDB files) and detect modellable elements -
i.e. protein, nucleic acid chains and ligands.

Classes
-------

Subdivision: Stored information on single type components (e.g. protein chain).
             Also provides basic functionality to align with a FASTA sequence
             where appropriate for modelling.

SourceStructure: Stores information on PDB structure and header and provides
                 functions to decompose into `Subdivisions`.
"""
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


class Subdivision:
    """
    Stores information on single type (i.e. protein, dna, rna, water ...)
    subset of residues in a structure (from a single PDB chain) and provides
    functions to perform basic alignment to FASTA sequences (where
    appropriate) and derive information for structural modelling.

    Parameters
    ----------
    residue_type: str
        What type of residues does the subdivision refer to - protein, dna,
        rna, water, ion or unknown.
    chain_id : str
        Which PDB chain is the subdivision derived from.
    structure: :class:`Structure <parmed.Structure>`
        Structure from which the subdivision is to be derived.
    """

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
        """
        Check if the `self.residue_type` indicates that the `Subdivision`
        represents a known polymer type - i.e. is one of 'protein, 'dna' or
        'rna'.

        Returns
        -------
        bool:
            Is this a protein or nucleic acid polymer?
        """

        if self.residue_type in ['protein', 'dna', 'rna']:
            return True

        return False

    def is_complete(self, strict=False):
        """
        Check if this `Subdivision` represents a complete polymer chain - i.e.
        has no gaps between bonded residues, or if they are present and not
        `strict` then do we have a sequence aligned to it that can be used to
        fill in gaps.

        Parameters
        ----------
        strict : bool
            Do we accept an aligned sequence as having no gaps?

        Returns
        -------
        bool:
            Is the sequence of `self.residues` gap free or has a sequence been
            aligned to allow gaps to be filled.

        """

        if not self.gaps:
            return True

        if not strict and self.aligned_sequence:
            return True

        return False

    def sequence(self, gap_char='-', seq_format='fasta'):
        """
        Return the sequence of the `self.residues` list of residues, with gaps
        between residues listed using `gap_char`. The output can be lists of
        one or three character residue codes or a FASTA format string (depending
        on the `seq_format` chosen).

        Parameters
        ----------
        gap_char : str
            Character to use when representing gaps in polymers (will be
            replicated three times in three letter code lists).
        seq_format : str
            Specifies the returned sequence format: 'letter' = single letter
            code list, 'rescode' = three letter code list and 'fasta' = FASTA
            format string.

        Returns
        -------
        list or str
            Either a list of residue codes (one or three letter) or a FASTA
            format string, depending on the `seq_format` chosen.

        """

        return residue_list_sequence(self.residues, seq_format=seq_format,
                                     gap_char=gap_char)

    def align_sequence(self, sequence):
        """
        Attempt to aligned `sequence` with `self.residues` using a basic string
        matching. If possible `self.aligned_sequence` is set to the subsection
        of `sequence` which overlaps the structural residues,
        `self.pre_sequence` any parts of sequence that come before the overlap
        and `self.pre_sequence` any that come after it.

        Parameters
        ----------
        sequence: str
            FASTA formatted string of the target sequence.

        Returns
        -------
        bool:
            Could the input sequence be aligned (using very simple string
            matching) with `self.residues` list of residues.
        """

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
        """
        Get the locations where residues would be added if a model was created
        to fill in gaps detected in the `self.residues` list and extra residues
        of the specified number were added pre and post structure.

        Parameters
        ----------
        pre_res : int
            Number of residues intended to be added before the structural
            residues.
        post_res : int
            Number of residues intended to be added after the structural
            residues.
        numbering : str
            Numbering scheme to use - default is to start at 1 (specified as
            'modeller' as that is what is used by the Modeller software)
            otherwise just use the residue number.

        Returns
        -------
        list:
            List of tuples containing first and last residue number given the
            specified numbering scheme.

        """

        if numbering == 'modeller':
            offset = 0
        else:
            offset = self.residues[0].number - 1

        sequence = '-'*pre_res + self.sequence() + '-'*post_res

        return [(m.start(0) + 1 + offset, m.end(0) + offset) for m
                in re.finditer('-+', sequence)]

    def target_sequence(self, pre_res=0, post_res=0):
        """
        Get FASTA format sequence form modelling - include specified regions
        before and after overlapping region.

        Parameters
        ----------
        pre_res : int
            Number of residues intended to be added before the structural
            residues.
        post_res : int
            Number of residues intended to be added after the structural
            residues.

        Returns
        -------
        str:
            FASTA style sequence of aligned sequence and selected length
            additions before and after the structure in `self.residues`

        """

        aligned_sequence = self.aligned_sequence

        if aligned_sequence:

            return (self.pre_sequence[:pre_res] + aligned_sequence +
                    self.post_sequence[:post_res])

        else:
            return ''

    def structure(self):
        """
        Copy the residues contained in the subdivision into a new
        `parmed.Structure`.

        Returns
        -------
        structure: :class:`Structure <parmed.Structure>`
            A `Structure` from `self.src_structure` which contains only
            the residues in this subversion.
        """

        struct = self.src_structure

        dataframe = struct.to_dataframe()
        residue_idxs = [res.idx for res in self.residues]

        return struct[dataframe.resid.isin(residue_idxs)]

    def modeller_selection(self, filename):
        """
        Create a modeller selection string for use in a PIR format alignment
        file to load the correct residues from a PDB created from the
        subdivision.

        Parameters
        ----------
        filename: str
            Name of the PDB file from which the coordinates will be read.

        Returns
        -------
        str:
            'structure' line for use in Modeller PIR alignment files.

        """

        if filename[-4:] == '.pdb':
            pdbname = filename[:-4]
        else:
            pdbname = filename

        chain_id = self.chain_id
        start = self.residues[0].resid
        end = self.residues[-1].resid

        return (f"structure:{pdbname}:{start}:{chain_id}:{end}:{chain_id}"
                f"::::")

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

        return (f"{self.__class__}( residue_idxs={first_idx}-{last_idx}, "
                f"residue_type={residue_type}, "
                f"aligned_sequence={aligned_sequence})")

    def __str__(self):

        return pprint.pformat(self.__dict__, indent=4)


class SourceStructure:
    """
    Stores a linked `parmed.structure` with header information and methods that
    produce a decomposition into single type subdivions + other functions
    designed to help model building.

    Parameters
    ----------
    structure: :class:`Structure <parmed.Structure>` or str
        Either a ParmEd structure or path to a structure file (only PDBs
        supported at present).
    header_filename :
        File from which to read a header (intended for when a structure is
        passed, if `structure` is a filename the header is read from that
        by default).

    """

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

        for residue in struct.residues:

            residue_type = get_residue_type(residue)

            if residue_type == 'unknown':
                self.nonstandard.append(residue)
                polymer_bond = check_residue_polymer_bonded(residue)
                if polymer_bond is not None:
                    residue_type = polymer_bond

            # TODO: Is there a less hacky way to do this?
            residue.residue_type = residue_type
            self.chain_map[residue.chain].append(residue)

        self.chains = self.chain_map.keys()
        self.decomposition = self.default_decomposition()

    def __repr__(self):

        return (f"{self.__class__}(structure={self.structure}, "
                f"source_file={self.source_file}, "
                f"header_file={self.header_filename})")

    def default_decomposition(self):
        """
        Create a dictionary of lists of `Subdivisions` describing contiguous
        regions of each chain which have the same residue type (protein, rna,
        dna, water, ion, unknown).

        Returns
        -------
        dict:
            keys = PDB chain ids and values = list of `Subdivisions`. Each
            `Subdivision` refers to a contiguous region of residues of the same
            type.

        """

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
                                                              directions=[1],
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

    def chain_sequence(self, chain_id, seq_format='fasta', src='structure'):
        """
        Return a representation of the residues contained in the selected PDB
        chain. The output can be lists of one or three character residue codes
        or a FASTA format string (depending on the `seq_format` chosen). All
        residues in the chain no matter what type are included.

        Parameters
        ----------
        chain_id : str
            PDB chain ID for which to produce sequence.
        seq_format : str
            Specifies the returned sequence format: 'letter' = single letter
            code list, 'rescode' = three letter code list and 'fasta' = FASTA
            format string.
        src

        Returns
        -------

        """

        if src not in ['structure', 'header']:
            raise RuntimeError("src must be 'structure', 'header'")

        decomposition = self.decomposition
        header_info = self.header

        if src == 'header':

            three_letter_seq = header_info.sequences[chain_id]

            if seq_format == 'fasta':

                sequence = convert_resname_list(three_letter_seq,
                                                seq_format=seq_format)

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
