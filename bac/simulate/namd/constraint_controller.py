from enum import Enum

from bac.utils.decorators import pathlike, advanced_property, positive_integer, positive_decimal, boolean, pdbcolumn

from bac.simulate.coding import Encodable


class BondType(Enum):
    all = 'all'
    water = 'water'
    none = 'none'


class ConstraintController(Encodable):
    def __init__(self, **kwargs):
        self.bond_constraint: BondConstraint = None
        self.harmonic_constraint: HarmonicConstraint = None
        self.atom_constraint: AtomConstraint = None
        self.extra_bonds = kwargs.get('extra_bonds')

    @pathlike
    def extra_bonds(self): pass


class BondConstraint(Encodable):
    def __init__(self, *, bonds=None, tolerance=None, iterations=None, die_on_error=None, use_settle=None):
        self.bonds = bonds
        self.tolerance = tolerance
        self.iterations = iterations
        self.die_on_error = die_on_error
        self.use_settle = use_settle

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


class HarmonicConstraint(Encodable):
    def __init__(self, *, exponent=None, scaling=None, **kwargs):
        self.exponent = exponent
        self.reference_position_file = kwargs.get('reference_position_file')
        self.force_constant_file = kwargs.get('force_constant_file')
        self.force_constant_column = kwargs.get('force_constant_column')
        self.scaling = scaling
        self.select_constraints = kwargs.get('select_constraints')
        self.select_constraints_x = kwargs.get('select_constraints_x')
        self.select_constraints_y = kwargs.get('select_constraints_y')
        self.select_constraints_z = kwargs.get('select_constraints_z')

    @positive_integer(default=0, validator=lambda _, x: x % 2 == 0)
    def exponent(self): pass

    @pathlike
    def reference_position_file(self): pass

    @pathlike
    def force_constant_file(self): pass

    @pdbcolumn
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


class AtomConstraint(Encodable):
    def __init__(self, *, calculate_fixed_atom_forces=None, **kwargs):
        self.calculate_fixed_atom_forces = calculate_fixed_atom_forces
        self.file = kwargs.get('file')
        self.column = kwargs.get('column')

    @boolean(default=False)
    def calculate_fixed_atom_forces(self): pass

    @pathlike(default=lambda s: s.run.coordinates)
    def file(self): pass

    @pdbcolumn
    def column(self): pass
