
from supproperty import (boolean, positive_decimal, float_vector,
                         positive_integer, decimal, pathlike)
from bac.utils.decorators import pdbcolumn
from bac.simulate.coding import Encodable


class PressureController(Encodable):
    def __init__(self, **kwargs):

        self.use_group_pressure = kwargs.get('use_group_pressure')
        self.use_flexible_cell = kwargs.get('use_flexible_cell')
        self.use_constant_ratio = kwargs.get('use_constant_ratio')
        self.use_constant_area = kwargs.get('use_constant_area')

    @boolean(default=False)
    def use_group_pressure(self): pass

    @boolean(default=False)
    def use_flexible_cell(self): pass

    @boolean(default=False)
    def use_constant_ratio(self): pass

    @boolean(default=False)
    def use_constant_area(self): pass


class BerendsenPressureCoupling(PressureController):
    def __init__(self, *, target=None, compressibility=None, relaxation_time=None, frequency=None, **kwargs):
        super(BerendsenPressureCoupling, self).__init__(**kwargs)

        self.target = target
        self.compressibility = compressibility
        self.relaxation_time = relaxation_time
        self.frequency = frequency

    @positive_decimal
    def target(self): pass

    @positive_decimal
    def compressibility(self): pass

    @positive_decimal
    def relaxation_time(self): pass

    @positive_integer(default=lambda self: self.simulation.non_bonded_controller.full_elect_frequency,
                      validator=lambda self, x: x % self.simulation.non_bonded_controller.full_elect_frequency == 0
                                             and x % self.simulation.non_bonded_controller.non_bonded_frequency == 0)
    def frequency(self): pass


class LangevinPistonPressureControl(PressureController):
    def __init__(self, *, target=None, period=None, decay=None, temperature=None, **kwargs):
        super(LangevinPistonPressureControl, self).__init__(**kwargs)

        self.target = target
        self.period = period
        self.decay = decay
        self.temperature = temperature
        self.surface_tension_target = kwargs.get('surface_tension_target')
        self.strain_rate = kwargs.get('strain_rate')
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

    @float_vector(default=(0, 0, 0), validator=lambda _, x: x.size == 3)
    def strain_rate(self): pass

    @boolean(default=False)
    def exclude_from_pressure(self): pass

    @pathlike
    def exclude_from_pressure_file(self): pass

    @pdbcolumn
    def exclude_from_pressure_column(self): pass
