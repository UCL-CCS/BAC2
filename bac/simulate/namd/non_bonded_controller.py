from enum import Enum

from bac.utils.decorators import *

from bac.simulate.coding import Encodable


class Interaction(Enum):
    one_two = '1-2'
    one_three = '1-3'
    one_four = '1-4'
    scaled_one_four = 'scaled1-4'
    none = 'none'


class MTSAlgorithm(Enum):
    impulse = 'impulse'
    verlet = 'verletI'
    constant = 'constant'
    naive = 'naive'


class LongSplitting(Enum):
    c1 = 'c1'
    c2 = 'c2'


class SplitPatchType(Enum):
    hydrogen = 'hydrogen'
    position = 'position'


class NonBondedController(Encodable):
    def __init__(self, **kwargs):

        self.cutoff = kwargs.get('cutoff')
        self.switching = kwargs.get('switching')
        self.switching_distance = kwargs.get('switching_distance')
        self.vdw_force_switching = kwargs.get('vdw_force_switching')
        self.excluded_interactions = kwargs.get('excluded_interactions')
        self.one_four_scaling = kwargs.get('one_four_scaling')
        self.dielectric_constant = kwargs.get('dielectric_constant')
        self.non_bonded_scaling = kwargs.get('non_bonded_scaling')
        self.vdw_geometric_sigma = kwargs.get('vdw_geometric_sigma')
        self.limit_distance = kwargs.get('limit_distance')
        self.LJ_correction = kwargs.get('LJ_correction')

        self.non_bonded_frequency = kwargs.get('non_bonded_frequency')
        self.full_elect_frequency = kwargs.get('full_elect_frequency')
        self.mts_algorithm = kwargs.get('mts_algorithm')
        self.long_splitting = kwargs.get('long_splitting')

        self.pairlist_distance = kwargs.get('pairlist_distance', self.cutoff)
        self.steps_per_cycle = kwargs.get('steps_per_cycle')
        self.split_patch = kwargs.get('split_patch')
        self.hgroup_cutoff = kwargs.get('hgroup_cutoff', 2.5)
        self.margin = kwargs.get('margin', 0.0)
        self.pairlist_min_processors = kwargs.get('pairlist_min_processors', 1)
        self.pairlists_per_cycle = kwargs.get('pairlists_per_cycle', 2)
        self.output_pairlists = kwargs.get('output_pairlists', 0)
        self.pairlist_shrink = kwargs.get('pairlist_shrink', 0.01)
        self.pairlist_grow = kwargs.get('pairlist_grow', 0.01)
        self.pairlist_trigger = kwargs.get('pairlist_trigger', 0.3)

        self.pme = None
        self.molly = None

        self.cell_basis_vector_1 = kwargs.get('cell_basis_vector_1')
        self.cell_basis_vector_2 = kwargs.get('cell_basis_vector_2')
        self.cell_basis_vector_3 = kwargs.get('cell_basis_vector_3')
        self.cell_origin = kwargs.get('cell_origin')

        self.extended_system: Path = kwargs.get('extended_system')
        self.xst_file = kwargs.get('xst_file')
        self.xst_frequency = kwargs.get('xst_frequency')

        self.wrap_water = kwargs.get('wrap_water')
        self.wrap_all = kwargs.get('wrap_all')
        self.wrap_nearest = kwargs.get('wrap_nearest')


    @positive_decimal
    def cutoff(self): pass

    @boolean(default=True)
    def switching(self): pass

    @positive_decimal(validator=lambda self, x: x <= self.cutoff)
    def switching_distance(self): pass

    @boolean(default=False)
    def vdw_force_switching(self): pass

    @advanced_property(type=Interaction)
    def excluded_interactions(self): pass

    @positive_decimal(default=1, validator=lambda _, x: x <= 1)
    def one_four_scaling(self): pass

    @decimal(default=1, validator=lambda _, x: x >= 1)
    def dielectric_constant(self): pass

    @positive_decimal(default=1)
    def non_bonded_scaling(self): pass

    @boolean(default=False)
    def vdw_geometric_sigma(self): pass

    @positive_decimal(default=0)
    def limit_distance(self): pass

    @boolean(default=False)
    def LJ_correction(self): pass

    @positive_integer(default=lambda self: self.non_bonded_frequency,
                      validator=lambda self, x: self.steps_per_cycle % x == 0, warn=True)
    def full_elect_frequency(self): pass

    @positive_integer(default=1, validator=lambda self, x: x % self.full_elect_frequency == 0)
    def non_bonded_frequency(self): pass

    @advanced_property(type=LongSplitting, default=LongSplitting.c1)
    def long_splitting(self): pass

    @advanced_property(type=MTSAlgorithm, default=MTSAlgorithm.impulse)
    def mts_algorithm(self): pass

    @positive_integer(default=20)
    def steps_per_cycle(self): pass

    @advanced_property(type=SplitPatchType, default=SplitPatchType.hydrogen)
    def split_patch(self): pass

    # Boundary Condition

    @advanced_property(type=tuple, default=(0, 0, 0))
    def cell_basis_vector_1(self): pass

    @advanced_property(type=tuple, default=(0, 0, 0))
    def cell_basis_vector_2(self): pass

    @advanced_property(type=tuple, default=(0, 0, 0))
    def cell_basis_vector_3(self): pass

    @advanced_property(type=tuple, default=(0, 0, 0))
    def cell_origin(self): pass

    @pathlike
    def extended_system(self): pass

    @extended_system.post_set_processing
    def extended_system(self, value):
        if value is not None:
            self.cell_basis_vector_1 = None
            self.cell_basis_vector_2 = None
            self.cell_basis_vector_3 = None
            self.cell_origin = None

    @pathlike(default=lambda self: self.simulation.output_name.with_suffix('.xst') if self.xst_frequency else None)
    def xst_file(self): pass

    @positive_integer
    def xst_frequency(self): pass

    @boolean(default=False)
    def wrap_water(self):pass

    @boolean(default=False)
    def wrap_all(self): pass

    @boolean(default=False)
    def wrap_nearest(self): pass


class PME(Encodable):
    def __init__(self, **kwargs):
        self.tolerance = kwargs.get('tolerance', 10e-6)
        self.interpolation_order = kwargs.get('interpolation_order', 4)
        self.grid_spacing = kwargs.get('grid_spacing', 1.0)
        self.grid_size = kwargs.get('grid_size')


class Molly(Encodable):
    def __init__(self, **kwargs):
        self.tolerance = kwargs.get('tolerance', 0.00001)
        self.iterations = kwargs.get('iterations', 100)