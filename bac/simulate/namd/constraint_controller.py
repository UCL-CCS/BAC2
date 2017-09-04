from enum import Enum

from bac.utils.decorators import *


class BondType(Enum):
    all = 'all'
    water = 'water'
    none = 'none'


class ConstraintController:
    def __init__(self, **kwargs):
        self.bond_constraint = None
        self.harmonic_constraint = None
        self.atom_constraint = None
        self.extra_bonds = kwargs.get('extra_bonds')

    @file
    def extra_bonds(self): pass


class BondConstraint:
    def __init__(self, **kwargs):
        self.bonds = kwargs.get('bonds')
        self.tolerance = kwargs.get('tolerance')
        self.iterations = kwargs.get('iterations')
        self.die_on_error = kwargs.get('die_on_error')
        self.use_settle = kwargs.get('use_settle')

    @advanced_property(type=BondType, default=BondType.none)
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
    def __init__(self, **kwargs):
        self.exponent = kwargs.get('exponent')
        self.reference_position_file = kwargs.get('reference_position_file')
        self.force_constant_file = kwargs.get('force_constant_file')
        self.force_constant_column = kwargs.get('force_constant_column')
        self.scaling = kwargs.get('scaling')
        self.select_constraints = kwargs.get('select_constraints')
        self.select_constraints_x = kwargs.get('select_constraints_x')
        self.select_constraints_y = kwargs.get('select_constraints_y')
        self.select_constraints_z = kwargs.get('select_constraints_z')

    @positive_integer(default=0, validator=lambda _, x: x % 2 == 0)
    def exponent(self): pass

    @file
    def reference_position_file(self): pass

    @file
    def force_constant_file(self): pass

    @column
    def force_constant_column(self): pass

    @positive_decimal(default=1.0)
    def scaling(self): pass

    @boolean(default=False)
    def select_constraints(self): pass

    @boolean(default=False)
    def select_constraints_x(self): pass

    @boolean(default=False)
    def select_constraints_y(self): pass

    @boolean(default=False)
    def select_constraints_z(self): pass


class AtomConstraint:
    def __init__(self, **kwargs):
        self.calculated_fixed_atom_forces = kwargs.get('calculated_fixed_atom_forces')
        self.file = kwargs.get('file')
        self.column = kwargs.get('column')

    @boolean(default=False)
    def calculated_fixed_atom_forces(self): pass

    @file(default=lambda s: s.run.coordinates)
    def file(self): pass

    @column
    def column(self): pass
