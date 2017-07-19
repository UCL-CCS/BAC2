from enum import Enum

from bac.utils.decorators import advanced_property, integer, decimal
from bac.simulate.gromacs.run_controller import Integrator


class TemperatureCouple(Enum):
    no = 'no'
    berendsen = 'berendsen'
    nose_hoover = 'nose-hoover'
    andersen = 'andersen'
    andersen_massive = 'andersen-massive'
    velocity_rescale = 'v-rescale'


class TemperatureController:

    def __init__(self, **kwargs):

        self.coupling_type = kwargs.get('coupling_type')
        self.coupling_frequency = kwargs.get('coupling_frequency')
        self.nose_hoover_chain_length = kwargs.get('nose_hoover_chain_length')
        self.coupling_time = kwargs.get('coupling_time')
        self.coupling_temperature = kwargs.get('coupling_temperature')

    @advanced_property(type=TemperatureCouple, default=TemperatureCouple.no)
    def coupling_type(self): pass

    def coupling_frequency_default(self):
        if self.run.integrator in (Integrator.md_vv,Integrator.md_vv_avek):
            return 1
        else:
            return self.run.non_bonded_controller.nstlist if self.run.non_bonded_controller.nstlist > 0 else 10

    @integer(default=coupling_frequency_default)
    def coupling_frequency(self): pass

    @integer(default=10, validator=lambda v, o: (v == 1) if (o.run.integrator is Integrator.md) else True)
    def nose_hoover_chain_length(self): pass

    @decimal
    def coupling_time(self): pass

    @decimal
    def coupling_temperature(self): pass


