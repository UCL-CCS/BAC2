from collections import OrderedDict
import numpy as np
from .assembly import Symmetry, BioTransform, BioUnit
from .sequence import convert_resname_list


class HeaderInfo(object):
    """

    Attributes
    ----------
    sequences : dict
        Sequence, as list of three letter residue codes, for each chain.
    unit_transforms : dict
        Biological unit transforms.
    symmetry :class:`Symmetry <bac.builder.utils.assembly.Symmetry>`
        Symmetry tensor.
    disulphides : list
        Disulphide bond descriptions.
    """

    def __init__(self):

        self.sequences = {}
        self.unit_transforms = {}
        self.symmetry = {}
        self.disulphides = {}

    def chain_sequence(self, chain_id, seq_format='letter'):
        """

        Parameters
        ----------
        chain_id: str
            Chain identifier for which the sequence is required
        seq_format: str
            Choice of format to return - 'letter' = list of single letter residue
            codes, 'fasta' = fasta style string of single letter codes

        Returns
        -------
        list or str
            Either list or string of single letter residue codes
        """

        sequences = self.sequences

        if chain_id not in sequences:
            if seq_format == 'letter'
                return []
            else:
                return ''

        sequence = sequences[chain_id]

        return convert_resname_list(sequence, seq_format=seq_format)

    def add_seq_res_line(self, line):
        """
        Parses the input SEQRES line to obtain a list of three letter residue
        codes which are added to the `self.sequences` dictionary with the key
        being the chain id (also obtained from the SEQRES information).

        Parameters
        ----------
        line : str
            SEQRES line from a PDB header
        """

        sequences = self.sequences

        if not line.startswith('SEQRES'):
            raise RuntimeError(f'Unable to interpret line as SEQRES:\n{line}')

        chain_id = line[11]
        residue_names = line[19:70].split()

        if chain_id in sequences:
            sequences[chain_id] += residue_names
        else:
            sequences[chain_id] = residue_names


def _parse_biomt(lines):
    """ Parse BIOMT lines and return an ordered dictionary of `BioUnits`

    Parameters
    ----------
    lines : list
        REMARK 350 lines from a PDB.

    Returns
    -------
    biological_units : dict
        Dictionary containing BioUnits derived from PDB header lines (keys are
        identifiers given in the source).
    """

    biological_units = OrderedDict()

    biomol_no = None

    for line in lines:

        if 'BIOMOLECULE' in line:
            # Start of a new biological unit
            # Record any previous & reset variables to be read
            # Example of line format:
            # REMARK 350 BIOMOLECULE: 1

            if biomol_no is not None:

                biological_units[biomol_no] = BioUnit(chains=chains,
                                                      source=origin,
                                                      nmer=nmer,
                                                      transforms=matrices)

            matrices = []
            current_matrix = []
            nmer = ''
            chains=[]

            biomol_no = int(line.split()[3])

        elif 'DETERMINED' in line:
            # Is this an author or software determined unit?
            # What order of multimer is indicated?
            # Example of line format:
            # REMARK 350 SOFTWARE DETERMINED QUATERNARY STRUCTURE: DIMERIC
            nmer = line.split(':')[1]
            origin = line.split()[2]

        elif 'APPLY THE FOLLOWING' in line:
            # Which chains does this apply to?
            # Example of line format:
            # REMARK 350 APPLY THE FOLLOWING TO CHAINS: A,B

            chain_txt = line.split(':')[1]
            chains = [x.strip() for x in chain_txt.split(',')]

        elif '350   BIOMT' in line:
            # Parse the transformation matrices, provided for each
            # additional copy of the system elements to be created
            # Matrix is aplit across three lines (BIOMT1, BIOMT2, BIOMT3)
            # Example of line format:
            # REMARK 350   BIOMT1   1  1.000000  0.000000  0.000000        0.00000

            columns = line.split()
            row_label = columns[2][-1]
            current_matrix.append(columns[4:8])

            # BioUnit requires a 3*4 matrix containing translation and
            # rotation information
            if row_label == '3':
                matrices.append(BioTransform(current_matrix))
                current_matrix = []

    biological_units[biomol_no] = BioUnit(chains=chains, transforms=matrices,
                                          source=origin, nmer=nmer)

    return biological_units


def _parse_symmetry(lines):
    """
    Parse symmetry information from REMARK 290 PDB header lines.

    Parameters
    ----------
    lines : str
        REMARK 290 lines from a PDB.

    Returns
    -------
    :class:`Symmetry <bac.builder.utils.assembly.Symmetry>`
        Tensor describing crystal symmetry
    """

    data = []

    for line in lines:
        if line.strip().startswith('REMARK 290   SMTRY'):
            data.append(line.split()[4:])

    tensor = np.asarray(data, dtype='f8')

    return Symmetry(tensor)


def read_pdb_header(pdb_filename):
    """
    Read a PDB file in order to obtain header information.

    Parameters
    ----------
    pdb_filename: str
        Path to PDB file from which to read header

    Returns
    -------
    :class:`HeaderInfo <bac.builder.utils.header.HeaderInfo>`
        Information on PDB structure obtained from its header.
    """

    header_info = HeaderInfo()

    biomt_lines = []
    symmetry_lines = []

    with open(pdb_filename, 'r') as pdb:

        for line in pdb:

            if line[:6] in ['ATOM  ', 'HETATM']:
                break

            elif line.startswith('REMARK 350'):
                biomt_lines.append(line)

            elif line.startswith('REMARK 290   SMTRY'):
                symmetry_lines.append(line)

            elif line.startswith('SEQRES'):
                header_info.add_seq_res_line(line)

    header_info.unit_transforms = _parse_biomt(biomt_lines)
    header_info.symmetry = _parse_symmetry(symmetry_lines)

    return header_info
