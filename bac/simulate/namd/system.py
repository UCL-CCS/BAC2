from enum import Enum

import parmed as pmd

from bac.utils.decorators import *
from bac.simulate.coding import Encodable


class Interaction(Enum):
    one_two = '1-2'
    one_three = '1-3'
    one_four = '1-4'
    scaled_one_four = 'scaled1-4'
    none = 'none'


class SplitPatchType(Enum):
    hydrogen = 'hydrogen'
    position = 'position'


class WaterModel(Enum):
    tip3p = 'tip3'
    tip4p = 'tip4'
    swm4_ndp = 'swm4'


class System(Encodable):
    def __init__(self, **kwargs):

        # Main system descriptors like topology, coordinates etc.

        self.topology = kwargs.get('topology')

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

        self.pme = kwargs.get('pme')
        self.molly = None

        # self.cell_basis_vector_1 = kwargs.get('cell_basis_vector_1')
        # self.cell_basis_vector_2 = kwargs.get('cell_basis_vector_2')
        # self.cell_basis_vector_3 = kwargs.get('cell_basis_vector_3')
        # self.cell_origin = kwargs.get('cell_origin')

        self.extended_system: Path = kwargs.get('extended_system')

        self.xst_file = kwargs.get('xst_file')
        self.xst_frequency = kwargs.get('xst_frequency')

        self.wrap_water = kwargs.get('wrap_water')
        self.wrap_all = kwargs.get('wrap_all')
        self.wrap_nearest = kwargs.get('wrap_nearest')

        # self.water_model = kwargs.get('water_model')

    _PRMTOP_NAME = 'complex.prmtop'

    @advanced_property(type=pmd.Structure)
    def topology(self): pass

    @topology.post_set_processing
    def topology(self):
        if self.topology.box_vectors is not None:
            print('Setting cell vectors.')
            bv = self.topology.box_vectors.value_in_unit(pmd.unit.angstrom)
            self.cell_basis_vector_1 = bv[0]
            self.cell_basis_vector_2 = bv[1]
            self.cell_basis_vector_3 = bv[2]

            water = next(res for res in self.topology.residues if res.name == 'WAT')
            if len(water) == 3:
                self.water_model = WaterModel.tip3p
            elif len(water) == 4:
                self.water_model = WaterModel.tip4p

            print('Setting', self.water_model)

    @pathlike(default=_PRMTOP_NAME)
    def parameters(self): pass

    @positive_decimal
    def cutoff(self): pass

    @boolean(default=lambda self: bool(self.switching_distance))
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

    @positive_integer(default=20)
    def steps_per_cycle(self): pass

    @advanced_property(type=SplitPatchType, default=SplitPatchType.hydrogen)
    def split_patch(self): pass

    # Boundary Condition

    @float_vector(default=(0, 0, 0), validator=lambda _, x: x.size == 3)
    def cell_basis_vector_1(self): pass

    @float_vector(default=(0, 0, 0), validator=lambda _, x: x.size == 3)
    def cell_basis_vector_2(self): pass

    @float_vector(default=(0, 0, 0), validator=lambda _, x: x.size == 3)
    def cell_basis_vector_3(self): pass

    @float_vector(default=(0, 0, 0), validator=lambda _, x: x.size == 3)
    def cell_origin(self): pass

    @pathlike
    def extended_system(self): pass

    @extended_system.post_set_processing
    def extended_system(self):
        if self.extended_system is not None:
            print('Deleting cell vectors.')
            self.cell_basis_vector_1 = None
            self.cell_basis_vector_2 = None
            self.cell_basis_vector_3 = None
            self.cell_origin = None

    @pathlike(default=lambda self: self.simulation.output_name.with_suffix('.xst') if self.xst_frequency else None)
    def xst_file(self): pass

    @positive_integer
    def xst_frequency(self): pass

    @boolean(default=False)
    def wrap_water(self): pass

    @boolean(default=False)
    def wrap_all(self): pass

    @boolean(default=False)
    def wrap_nearest(self): pass

    @advanced_property(type=WaterModel, default=WaterModel.tip3p)
    def water_model(self): pass

    def encode(self, path=None, suffix=None):

        if path is not None:
            self.topology.save(self._PRMTOP_NAME)

        return super().encode(path=path, suffix=suffix)


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
