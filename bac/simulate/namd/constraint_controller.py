from enum import Enum

from bac.utils.decorators import *


class ConstraintController:
    def __init__(self):
        self.bond_constraint = None
        self.harmonic_constraint = None
        self.atom_constraint = None
        self.extra_bonds = None


class BondType(Enum):
    all = 'all'
    water = 'water'


class BondConstraint:
    def __init__(self, **kwargs):
        self.bonds = kwargs.get('bonds')
        self.tolerance = kwargs.get('tolerance')
        self.iterations = kwargs.get('iterations')
        self.die_on_error = kwargs.get('die_on_error')
        self.use_settle = kwargs.get('use_settle')

    @advanced_property(type=(BondType, str))
    def bonds(self): pass

    @positive_decimal(default=1.0e-8)
    def tolerance(self): pass

    @positive_integer
    def iterations(self): pass

    @boolean(default=True)
    def die_on_error(self): pass

    @boolean(default=True)
    def use_settle(self): pass


class HarmonicConstraint:
    # TODO: implement
    pass


class AtomConstraint:
    def __init__(self, **kwargs):
        self.calculated_fixed_atom_forces = kwargs.get('calculated_fixed_atom_forces', False)
        self.file = kwargs.get('file')
        self.column = kwargs.get('column', 'O')


class ExtraBonds:
    def __init__(self, file=None):
        self.file = file
