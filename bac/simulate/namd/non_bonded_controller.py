from enum import Enum

from bac.utils.decorators import *


class Interaction(Enum):
    one_two = '1-2'
    one_three = '1-3'
    one_four = '1-4'
    scaled_one_four = 'scaled1-4'
    none = 'none'


class MTSAlgorithm(Enum):
    impulse_verlet = 'impulse/verletI'
    constant_naive = 'constant/naive'


class LongSplitting(Enum):
    c1 = 'c1'
    c2 = 'c2'


class SplitPatchType(Enum):
    hydrogen = 'hydrogen'
    position = 'position'


class NonBondedController:

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

    @positive_decimal
    def cutoff(self): pass

    @boolean(default=True)
    def switching(self): pass

    @positive_decimal(validator=lambda x, s: x <= s.cutoff)
    def switching_distance(self): pass

    @boolean(default=False)
    def vdw_force_switching(self): pass

    @advanced_property(type=Interaction)
    def excluded_interactions(self): pass

    @positive_decimal(default=1, validator=lambda x, s: x <= 1)
    def one_four_scaling(self): pass

    @decimal(default=1, validator=lambda x: x >= 1)
    def dielectric_constant(self): pass

    @positive_decimal(default=1)
    def non_bonded_scaling(self): pass

    @boolean(default=False)
    def vdw_geometric_sigma(self): pass

    @positive_decimal(default=0)
    def limit_distance(self): pass

    @boolean(default=False)
    def LJ_correction(self): pass


    @positive_integer(default=lambda x: x.non_bonded_frequency,
                      validator=lambda x, s: s.steps_per_cycle % x == 0, warn=True)
    def full_elect_frequency(self): pass

    @positive_integer(default=1, validator=lambda x, s: x % s.full_elect_frequency == 0)
    def non_bonded_frequency(self): pass

    @advanced_property(type=LongSplitting, default=LongSplitting.c1)
    def long_splitting(self): pass

    @advanced_property(type=MTSAlgorithm, default=MTSAlgorithm.impulse_verlet)
    def mts_algorithm(self): pass




    @positive_integer(default=20)
    def steps_per_cycle(self): pass

    @advanced_property(type=SplitPatchType, default=SplitPatchType.hydrogen)
    def split_patch(self): pass


class PME:
    def __init__(self, **kwargs):
        self.tolerance = kwargs.get('tolerance', 10e-6)
        self.interpolation_order = kwargs.get('interpolation_order', 4)
        self.grid_spacing = kwargs.get('grid_spacing', 1.0)
        self.grid_size = kwargs.get('grid_size')


class Molly:
    def __init__(self, **kwargs):
        self.tolerance = kwargs.get('tolerance', 0.00001)
        self.iterations = kwargs.get('iterations', 100)


