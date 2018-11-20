import numpy as np


class Symmetry(object):
    """Symmetry operations

    Parameters
    ----------
    tensor : np.array, shape=(3*n_symmetry, 4)
    """
    def __init__(self, tensor):
        self.data = np.asarray(tensor, dtype='f8')
        if self.data.ndim != 2 or self.data.shape[1] != 4:
            raise ValueError('Tensor shape must be (3*n_symmetry, 4)')


class BioTransform(object):
    """Biological unit transform

    Parameters
    ----------
    matrix: np.array, shape=(3, 4)
        Matrix describing rotation and translation as described in PDB headers
    """
    def __init__(self, matrix):
        self.data = np.asarray(matrix, np.float)
        if ((self.data.ndim != 2)
           or (self.data.shape[0] != 3)
           or (self.data.shape[1] != 4)):

            raise ValueError('Matrix must be (3, 4)')

    def transform_coords(self, coords):
        """
        Apply the transform described by the `self.data` matrix to the provided
        coordinates.

        Parameters
        ----------
        coords: np.array
            Atomic coordinates to be transformed.

        Returns
        -------
        np.array
            Coordinates after the application of the transform.
        """

        coords4 = np.ones((len(coords), 4), dtype=np.float)
        coords4[:, 0:3] = coords

        return coords4.dot(self.data.T)


class BioUnit(object):
    """Biological unit information

    Parameters
    ----------
    transforms : list
        Transforms to apply to chains in the provided coordinates to create the
        biological unit.
    chains : list
        Chains to which the transforms applies.
    source : str
        Description of the source of the biological unit determination (author
        or software).
    nmer   : str
        Description of the target multimer.
    """

    def __init__(self, transforms=[], chains=[],
                 source="", nmer="",  **kwargs):

        self.source = source
        self.chains = chains
        self.nmer = nmer
        self.transforms = transforms
