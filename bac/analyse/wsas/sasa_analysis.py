import numpy as np


def atom_contribution_nm(atom_type, sas, params):
    """
    Calculate the contribution for a given atom to the estimation of the
    configurational entropy that would be computed via normal mode analysis.
    Formula from:
        Wang, J., & Hou, T. (2012). Develop and Test a Solvent Accessible
        Surface Area-Based Model in Conformational Entropy Calculations.
        Journal of Chemical Information and Modeling,
        52(5), 1199–1212. http://doi.org/10.1021/ci300064d

    contribution = weight * (sas + k * bsa)
    bsa = (4 * pi * (radius + probe_radius)^2) - sas

    sas = solvent accessible surface
    bsa = buried surface area

    Note: Units are cal/mol (NOT kcal/mol)

    Parameters
    ----------
    atom_type : str
        Atom type of selected atom.
    sas: float
        Solvent accessible surface area for selected atom.
    params : dict
        Dictionary containing sas to normal mode configurational entropy
        estimate parameters.

    Returns
    -------
    float
        Atomic contribution to the configurational entropy

    """

    probe_radius = params['rprobe']
    atom_type_info = params['params'][atom_type]
    k = params['k']

    total_atom_surf = 4 * np.pi * (atom_type_info['radius'] + probe_radius) ** 2
    bsa = total_atom_surf - sas

    return atom_type_info['weight'] * (sas + k * bsa)

def nm_component_calc(atom_values, temperature, intercept):

    #Factor of 1000 converts cal/mol to kcal/mol
    component_nm = temperature/1000 * (atom_values.sum() + intercept)

    return component_nm