from bac.utils.decorators import (boolean, positive_decimal, pathlike, pdbcolumn,
                                  positive_integer, decimal, advanced_property)

from bac.simulate.coding import Encodable
from bac.utils.pdb import PDBColumn
import bac.simulate.base_thermostat as base


class TemperatureController:
    pass


class LangevinDynamics(TemperatureController, Encodable):
    def __init__(self, *, temperature, damping=None, applies_to_hydrogen=None, file=None, column=None):
        self.temperature = temperature
        self.damping = damping
        self.applies_to_hydrogen = applies_to_hydrogen
        self.file = file
        self.column = column

    @positive_decimal
    def temperature(self): pass

    @positive_decimal
    def damping(self): pass

    @boolean(default=lambda self: True if self.damping is not None else None)
    def applies_to_hydrogen(self): pass

    @pathlike(default=lambda self: self.simulation.coordinates if self.damping is None else None)
    def file(self): pass

    @advanced_property(type=PDBColumn, default=lambda self: PDBColumn.O if self.damping is None else None)
    def column(self): pass


class TemperatureCoupling(TemperatureController, Encodable):

    def __init__(self, *, temperature, file=None, column=None):
        self.temperature = temperature
        self.file = file
        self.column = column

    @positive_decimal
    def temperature(self): pass

    @pathlike(default=lambda self: self.simulation.coordinates)
    def file(self): pass

    @pdbcolumn
    def column(self): pass


class VelocityRescale(base.VelocityRescale):
    def __init__(self, *, temperature, frequency):
        super(VelocityRescale, self).__init__(temperature=temperature, frequency=frequency)

    @positive_decimal
    def temperature(self): pass


class VelocityRescaling(TemperatureController, Encodable):
    def __init__(self, *, temperature, frequency):
        self.temperature = temperature
        self.frequency = frequency

    @positive_decimal
    def temperature(self): pass

    @positive_integer
    def frequency(self): pass

    def __bool__(self):
        return self.frequency is not None


class VelocityReassignment(TemperatureController, Encodable):
    def __init__(self, *, temperature=None, frequency, increment=None, hold_at=None):
        self.temperature = temperature
        self.frequency = frequency
        self.increment = increment
        self.hold_at = hold_at

    @positive_decimal(default=lambda self: self.simulation.temperature)
    def temperature(self): pass

    @positive_integer
    def frequency(self): pass

    @decimal(default=0)
    def increment(self): pass

    @positive_decimal
    def hold_at(self): pass

    def __bool__(self):
        return self.frequency is not None


class LoweAndersenDynamics(TemperatureController, Encodable):
    def __init__(self, *, temperature, cutoff=None, rate=None):
        self.temperature = temperature
        self.cutoff = cutoff
        self.rate = rate

    @positive_decimal
    def temperature(self): pass

    @positive_decimal(default=2.7)
    def cutoff(self): pass

    @positive_decimal(default=50)
    def rate(self): pass
