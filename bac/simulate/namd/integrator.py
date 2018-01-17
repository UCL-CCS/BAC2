from enum import Enum

from bac.utils.decorators import advanced_property, positive_integer, positive_decimal, boolean, pathlike, pdbcolumn
from bac.simulate.coding import Encodable


class Integrator:
    def __init__(self, timestep=None, first_step=None):
        self.timestep = timestep
        self.first_step = first_step

    @positive_decimal(default=1.0)
    def timestep(self): pass

    @positive_integer(default=0)
    def first_step(self): pass


class VerletIntegrator(Integrator, Encodable):
    def __init__(self, *, timestep=None, first_step=None):
        super(VerletIntegrator, self).__init__(timestep, first_step)


class LongSplitting(Enum):
    c1 = 'c1'
    c2 = 'c2'


class MTSAlgorithm(Enum):
    impulse = 'impulse'
    constant = 'constant'


class MTSIntegrator(Integrator, Encodable):
    def __init__(self, timestep=None, first_step=None, non_bonded_frequency=None,
                 electrostatic_frequency=None, long_splitting=None, algorithm=None):
        super(MTSIntegrator, self).__init__(timestep, first_step)
        self.timestep = timestep
        self.non_bonded_frequency = non_bonded_frequency
        self.electrostatic_frequency = electrostatic_frequency
        self.long_splitting = long_splitting
        self.algorithm = algorithm

    @positive_integer(default=lambda self: self.non_bonded_frequency,
                      validator=lambda self, x: self.simulation.non_bonded_controller.steps_per_cycle % x == 0, warn=True)
    def electrostatic_frequency(self): pass

    @positive_integer(default=1, validator=lambda self, x: x % self.electrostatic_frequency == 0)
    def non_bonded_frequency(self): pass

    @advanced_property(type=LongSplitting, default=LongSplitting.c1)
    def long_splitting(self): pass

    @advanced_property(type=MTSAlgorithm, default=MTSAlgorithm.impulse)
    def algorithm(self): pass


class LangevinIntegrator(MTSIntegrator, Encodable):
    def __init__(self, timestep=None, first_step=None, temperature=None, damping=None, applies_to_hydrogen=None, file=None, column=None, **kwargs):
        super(LangevinIntegrator, self).__init__(timestep, first_step, **kwargs)
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

    @pdbcolumn(default=lambda self: 'O' if self.damping is None else None)
    def column(self): pass
