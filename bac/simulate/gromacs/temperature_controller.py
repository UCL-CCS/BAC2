from enum import Enum
from functools import partial

import numpy as np

from bac.utils.decorators import advanced_property, integer
from bac.simulate.gromacs.integrator import Integrator
from bac.simulate.coding import Encodable


class TemperatureCouplingType(Enum):
    no = 'no'
    berendsen = 'berendsen'
    nose_hoover = 'nose-hoover'
    andersen = 'andersen'
    andersen_massive = 'andersen-massive'
    velocity_rescale = 'v-rescale'

    @classmethod
    def _missing_(cls, value):
        if value is False:
            return cls.no
        else:
            super()._missing_(value)


class TemperatureController(Encodable):
    def __init__(self, **kwargs):

        self.type = kwargs.get('type')
        self.frequency = kwargs.get('frequency')
        self.nose_hoover_chain_length = kwargs.get('nose_hoover_chain_length')
        self.groups = kwargs.get('groups')
        self.time = kwargs.get('time')
        self.temperature = kwargs.get('temperature')

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return (self.type == other.type and
                    self.frequency == other.frequency and
                    self.nose_hoover_chain_length == other.nose_hoover_chain_length and
                    self.time == other.time and
                    self.temperature == other.temperature)

        return False

    @advanced_property(type=TemperatureCouplingType, default=TemperatureCouplingType.no)
    def type(self): pass

    def _frequency_default(self):
        if self.simulation.integrator in (Integrator.md_vv, Integrator.md_vv_avek):
            return 1
        else:
            return self.simulation.non_bonded_controller.neighbor_list_update_frequency \
                if self.simulation.non_bonded_controller.neighbor_list_update_frequency > 0 else 10

    @integer(default=_frequency_default)
    def frequency(self): pass

    @integer(default=10, validator=lambda self, v: (v == 1) if (self.simulation.integrator is Integrator.md) else True)
    def nose_hoover_chain_length(self): pass

    @advanced_property(type=list, default=[])
    def groups(self): pass

    @advanced_property(type=partial(np.array, dtype=np.float), default=[],
                       validator=lambda self, x: self.groups.size == x.size)
    def time(self): pass

    @advanced_property(type=partial(np.array, dtype=np.float), default=[],
                       validator=lambda self, x: self.groups.size == x.size)
    def temperature(self): pass
