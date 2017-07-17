from enum import Enum

from bac.utils.decorators import *


class LongSplitting(Enum):
    c1 = 'c1'
    c2 = 'c2'


class NonBondedController:

    def __init__(self, **kwargs):

        self.cutoff = kwargs.get('cutoff')
        self.switching = kwargs.get('switching', True)
        self.switching_distance = kwargs.get('switching_distance')
        self.vdw_force_switching = kwargs.get('vdw_force_switching', False)
        self.excluded_interactions = kwargs.get('excluded_interactions')
        self.one_four_scaling = kwargs.get('one_four_scaling', 1.0)
        self.dielectric_constant = kwargs.get('dielectric_constant', 1.0)
        self.non_bonded_scaling = kwargs.get('non_bonded_scaling', 1.0)
        self.vdw_geometric_sigma = kwargs.get('vdw_geometric_sigma', False)
        self.limit_distance = kwargs.get('limit_distance', 0.0)
        self.LJ_correction = kwargs.get('LJ_correction', False)

        self.non_bonded_frequency = kwargs.get('non_bonded_frequency')
        self.full_elect_frequency = kwargs.get('full_elect_frequency')
        self.mts_algorithm = kwargs.get('mts_algorithm', 'impulse/verlet')
        self.long_splitting = kwargs.get('long_splitting', 'c1')

        self.pairlist_distance = kwargs.get('pairlist_distance', self.cutoff)
        self.steps_per_cycle = kwargs.get('steps_per_cycle')
        self.split_patch = kwargs.get('split_patch', 'hydrogen')
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

    @advanced_property(type=(LongSplitting, str), default="c1")
    def long_splitting(self): pass

    @positive_integer(default=20)
    def steps_per_cycle(self): pass

    @positive_integer(default=lambda x: x.non_bonded_frequency,
                      validator=lambda x, s: s.steps_per_cycle % x == 0, warn=True)
    def full_elect_frequency(self): pass

    @positive_integer(default=1, validator=lambda x, s: x % s.full_elect_frequency == 0)
    def non_bonded_frequency(self): pass


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


