from enum import Enum

from bac.utils.decorators import advanced_property, decimal, positive_integer, positive_decimal, boolean


class CutoffSchemeType(Enum):
    verlet = 'Verlet'
    group = 'group'


class NeighborListType(Enum):
    grid = 'grid'
    simple = 'simple'


class PeriodicBoundaryConditionType(Enum):
    xyz = 'xyz'
    xy = 'xy'
    no = 'no'


class CoulombType(Enum):
    cutoff = 'Cut-off'
    ewald = 'Ewald'
    pme = 'PME'
    p3m_ad = 'P3M-AD'
    reaction_field = 'Reaction-Field'
    generalized_reaction_field = 'Generalized-Reaction-Field'
    reaction_field_zero = 'Reaction-Field-zero'
    shift = 'Shift'
    encad_shift = 'Encad-Shift'
    switch = 'Switch'
    user = 'User'
    pme_switch = 'PME-Switch'
    pme_user = 'PME-User'
    pme_user_switch = 'PME-User-Switch'


class CoulombModifierType(Enum):
    potential_shift_verlet = 'Potential-shift-Verlet'
    potential_shift = 'Potential-shift'
    none = 'None'


class VanDerWaalsType(Enum):
    cutoff = 'Cut-off'
    pme = 'PME'
    shift = 'Shift'
    switch = 'Switch'
    encad_shift = 'Encad-Shift'
    user = 'user'


class VanDerWaalsModifierType(Enum):
    potential_shift_verlet = 'Potential-shift-Verlet'
    potential_shift = 'Potential-shift'
    none = 'None'
    force_switch = 'Force-switch'
    potential_switch = 'Potential-switch'


class LongRangeDispersionCorrectionType(Enum):
    no = 'no'
    energy = 'Ener'
    energy_and_pressure = 'EnerPres'


class NonBondedController:

    def __init__(self, **kwargs):
        self.cutoff = kwargs.get('cutoff')
        self.cutoff_scheme = kwargs.get('cutoff_scheme')
        self.neighbor_list_update_frequency = kwargs.get('neighbour_list_update_frequency')
        self.neighbor_list_type = kwargs.get('neighbor_list_type')

        self.periodic_boundary_condition = kwargs.get('periodic_boundary_condition')
        self.periodic_molecules = kwargs.get('periodic_molecules')

        self.verlet_buffer_tolerance = kwargs.get('verlet_buffer_tolerance')

        self.coulomb_type = kwargs.get('coulomb_type')
        self.coulomb_modifier = kwargs.get('coulomb_modifier')
        self.coulomb_switch_cutoff = kwargs.get('coulomb_switch_cutoff')
        self.coulomb_cutoff = kwargs.get('coulomb_cutoff')
        self.dielectric = kwargs.get('dielectric')
        self.relative_dielectric = kwargs.get('relative_dielectric')

        self.van_der_waals_type = kwargs.get('van_der_waals_type')
        self.van_der_waals_modifier = kwargs.get('van_der_waals_modifier')
        self.van_der_waals_switch_cutoff = kwargs.get('van_der_waals_switch_cutoff')
        self.van_der_waals_cutoff = kwargs.get('van_der_waals_cutoff')
        self.correct_long_range_dispersion = kwargs.get('correct_long_range_dispersion')

    @positive_decimal(default=1.0)
    def cutoff(self): pass

    @advanced_property(type=CutoffSchemeType, default=CutoffSchemeType.verlet)
    def cutoff_scheme(self): pass

    @positive_integer(default=10)
    def neighbor_list_update_frequency(self): pass

    @advanced_property(type=NeighborListType, default=NeighborListType.grid)
    def neighbor_list_type(self): pass

    @advanced_property(type=PeriodicBoundaryConditionType, default=PeriodicBoundaryConditionType.xyz)
    def periodic_boundary_condition(self): pass

    @boolean(default=False)
    def periodic_molecules(self): pass

    @positive_decimal(default=0.005)
    def verlet_buffer_tolerance(self): pass

    @advanced_property(type=CoulombType, default=CoulombType.cutoff)
    def coulomb_type(self): pass

    @advanced_property(type=CoulombModifierType, default=CoulombModifierType.potential_shift_verlet)
    def coulomb_modifier(self): pass

    @decimal(default=0)
    def coulomb_switch_cutoff(self): pass

    @decimal(default=1)
    def coulomb_cutoff(self): pass

    @decimal(default=1)
    def dielectric(self): pass

    @decimal(default=0)
    def relative_dielectric(self): pass

    @advanced_property(type=VanDerWaalsType, default=VanDerWaalsType.cutoff)
    def van_der_waals_type(self): pass

    @advanced_property(type=VanDerWaalsModifierType, default=VanDerWaalsModifierType.potential_shift_verlet)
    def van_der_waals_modifier(self): pass

    @decimal(default=0)
    def van_der_waals_switch_cutoff(self): pass

    @decimal(default=1)
    def van_der_waals_cutoff(self): pass

    @advanced_property(type=LongRangeDispersionCorrectionType, default=LongRangeDispersionCorrectionType.no)
    def correct_long_range_dispersion(self): pass
