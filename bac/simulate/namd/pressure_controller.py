from bac.utils.pdb import PDBColumn
from bac.utils.decorators import (boolean, positive_decimal,
                                  positive_integer, decimal, file, column)


class PressureController:
    def __init__(self, **kwargs):

        self.use_group_pressure = kwargs.get('use_group_pressure')
        self.use_flexible_cell = kwargs.get('use_flexible_cell')
        self.use_constant_ratio = kwargs.get('use_constant_ratio')
        self.use_constant_area = kwargs.get('use_constant_area')
        self.berendsen = None
        self.langevin = None

    @boolean(default=False)
    def use_group_pressure(self): pass

    @boolean(default=False)
    def use_flexible_cell(self): pass

    @boolean(default=False)
    def use_constant_ratio(self): pass

    @boolean(default=False)
    def use_constant_area(self): pass


class BerendsenPressureCoupling:
    def __init__(self, **kwargs):
        self.target = kwargs.get('target')
        self.compressibility = kwargs.get('compressibility')
        self.relaxation_time = kwargs.get('relaxation_time')
        self.frequency = kwargs.get('frequency')

    @positive_decimal
    def target(self): pass

    @positive_decimal
    def compressibility(self): pass

    @positive_decimal
    def relaxation_time(self): pass

    @positive_integer(default=lambda self: self.run.non_bonded_controller.full_elect_frequency,
                      validator=lambda self, x: x % self.run.non_bonded_controller.full_elect_frequency == 0
                                             and x % self.run.non_bonded_controller.non_bonded_frequency == 0)
    def frequency(self): pass


class LangevinPistonPressureControl:
    def __init__(self, **kwargs):
        self.target = kwargs.get('target')
        self.period = kwargs.get('period')
        self.decay = kwargs.get('decay')
        self.temperature = kwargs.get('temperature')
        self.surface_tension_target = kwargs.get('surface_tension_target')
        self.strain_rate = kwargs.get('strain_rate', (0.0, 0.0, 0.0))
        self.exclude_from_pressure = kwargs.get('exclude_from_pressure')
        self.exclude_from_pressure_file = kwargs.get('exclude_from_pressure_file')
        self.exclude_from_pressure_column = kwargs.get('exclude_from_pressure_column')

    @positive_decimal
    def target(self): pass

    @positive_decimal
    def period(self): pass

    @positive_decimal
    def decay(self): pass

    @positive_decimal
    def temperature(self): pass

    @decimal(default=0.0)
    def surface_tension_target(self): pass

    @boolean(default=False)
    def exclude_from_pressure(self): pass

    @file
    def exclude_from_pressure_file(self): pass

    @column
    def exclude_from_pressure_column(self): pass
