from enum import Enum

from bac.utils.decorators import advanced_property, decimal, positive_integer, positive_decimal, boolean


class CutoffSchemeType(Enum):
    """Neighbor search cutoff scheme

    :param verlet: Generate a pair list with buffering. The buffer size is automatically set based on
    verlet-buffer-tolerance, unless this is set to -1, in which case rlist will be used. This option
    has an explicit, exact cut-off at rvdw equal to rcoulomb, unless PME or Ewald is used, in which
    case rcoulomb > rvdw is allowed. Currently only cut-off, reaction-field, PME or Ewald electrostatics
    and plain LJ are supported. Some gmx mdrun functionality is not yet supported with the Verlet scheme,
    but gmx grompp checks for this. Native GPU acceleration is only supported with Verlet. With GPU-accelerated
    PME or with separate PME ranks, gmx mdrun will automatically tune the CPU/GPU load balance by scaling rcoulomb
    and the grid spacing. This can be turned off with mdrun -notunepme. Verlet is faster than group when there is
    no water, or if group would use a pair-list buffer to conserve energy.
    :param group:
    """
    verlet = 'Verlet'
    group = 'group'


class NeighborListType(Enum):
    grid = 'grid'
    simple = 'simple'


class PeriodicBoundaryConditionType(Enum):
    xyz = 'xyz'
    xy = 'xy'
    no = 'no'

    @classmethod
    def _missing_(cls, value):
        if value is False:
            return cls.no
        else:
            super()._missing_(value)


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

    @classmethod
    def _missing_(cls, value):
        if value is None:
            return cls.none
        else:
            super()._missing_(value)


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

    @classmethod
    def _missing_(cls, value):
        if value is None:
            return cls.none
        else:
            super()._missing_(value)


class LongRangeDispersionCorrectionType(Enum):
    no = 'no'
    energy = 'Ener'
    energy_and_pressure = 'EnerPres'

    @classmethod
    def _missing_(cls, value):
        if value is False:
            return cls.no
        else:
            super()._missing_(value)


class EwaldCombinationType(Enum):
    geometric = 'Geometric'
    lorentz_berthelot = 'Lorentz-Berthelot'


class EwaldGeometryType(Enum):
    three_dimensions = '3d'
    three_dimensions_corrected = '3dc'


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

        self.fourier_spacing = kwargs.get('fourier_spacing')
        self.fourier_magnitude_x = kwargs.get('fourier_magnitude_x')
        self.fourier_magnitude_y = kwargs.get('fourier_magnitude_y')
        self.fourier_magnitude_z = kwargs.get('fourier_magnitude_z')
        self.pme_order = kwargs.get('pme_order')
        self.ewald_tolerance_coulomb = kwargs.get('ewald_tolerance_coulomb')
        self.ewald_tolerance_vdw = kwargs.get('ewald_tolerance_vdw')
        self.ewald_combination_rule = kwargs.get('ewald_combination_rule')
        self.ewald_geometry = kwargs.get('ewald_geometry')
        self.ewald_epsilon_surface = kwargs.get('ewald_epsilon_surface')

    @positive_decimal(default=1.0)
    def cutoff(self): pass

    @advanced_property(type=CutoffSchemeType, default=CutoffSchemeType.verlet,
                       validator=lambda _, x: x is CutoffSchemeType.verlet, warn=True,
                       warning_message='`CutoffSchemeType.group` was the only cut-off treatment scheme '
                                       'before version 4.6, and is deprecated in 5.0.')
    def cutoff_scheme(self):
        """ Cutoff scheme

        """
        pass

    @positive_integer(default=10)
    def neighbor_list_update_frequency(self): pass

    @advanced_property(type=NeighborListType, default=NeighborListType.grid)
    def neighbor_list_type(self): pass

    @advanced_property(type=(PeriodicBoundaryConditionType, bool, str), default=PeriodicBoundaryConditionType.xyz)
    def periodic_boundary_condition(self): pass

    @boolean(default=False)
    def periodic_molecules(self): pass

    @positive_decimal(default=0.005)
    def verlet_buffer_tolerance(self): pass

    @advanced_property(type=CoulombType, default=CoulombType.cutoff)
    def coulomb_type(self): pass

    @advanced_property(type=(CoulombModifierType, type(None), str), default=CoulombModifierType.potential_shift_verlet)
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

    @advanced_property(type=(VanDerWaalsModifierType, type(None), str),
                       default=VanDerWaalsModifierType.potential_shift_verlet)
    def van_der_waals_modifier(self): pass

    @decimal(default=0)
    def van_der_waals_switch_cutoff(self): pass

    @decimal(default=1)
    def van_der_waals_cutoff(self): pass

    @advanced_property(type=(LongRangeDispersionCorrectionType, bool, str), default=LongRangeDispersionCorrectionType.no)
    def correct_long_range_dispersion(self): pass

    @decimal(default=0.12)
    def fourier_spacing(self): pass

    @positive_integer(default=0)
    def fourier_magnitude_x(self): pass

    @positive_integer(default=0)
    def fourier_magnitude_y(self): pass

    @positive_integer(default=0)
    def fourier_magnitude_z(self): pass

    @positive_integer(default=4, validator=lambda _, x: x in [4, 6, 8, 10])
    def pme_order(self): pass

    @positive_decimal(default=1e-5)
    def ewald_tolerance_coulomb(self): pass

    @positive_decimal(default=1e-3)
    def ewald_tolerance_vdw(self): pass

    @advanced_property(type=EwaldCombinationType, default=EwaldCombinationType.geometric)
    def ewald_combination_rule(self): pass

    @advanced_property(type=EwaldGeometryType, default=EwaldGeometryType.three_dimensions)
    def ewald_geometry(self): pass

    @decimal(default=0)
    def ewald_epsilon_surface(self): pass
