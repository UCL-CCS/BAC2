from enum import Enum

from bac.utils.decorators import advanced_property, integer, decimal


class Integrator(Enum):
    md = 'md'
    md_vv = 'md-vv'
    md_vv_avek = 'md-vv-avek'
    sd = 'sd'
    bd = 'bd'
    steep = 'steep'
    cg = 'cg'
    l_bfgs = 'l-bfgs'
    nm = 'nm'
    tpi = 'tpi'
    tpic = 'tpic'

class TemperatureCouplingType(Enum):
    no = 'no'
    berendsen = 'berendsen'
    nose_hoover = 'nose-hoover'
    andersen = 'andersen'
    andersen_massive = 'andersen-massive'
    velocity_rescale = 'v-rescale'


class TemperatureController:
    def __init__(self, **kwargs):

        self.type = kwargs.get('type')
        self.frequency = kwargs.get('frequency')
        self.nose_hoover_chain_length = kwargs.get('nose_hoover_chain_length')
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

    def frequency_default(self):
        if self.run.integrator in (Integrator.md_vv,Integrator.md_vv_avek):
            return 1
        else:
            return self.run.non_bonded_controller.nstlist if self.run.non_bonded_controller.nstlist > 0 else 10

    @integer(default=frequency_default)
    def frequency(self): pass

    @integer(default=10, validator=lambda v, o: (v == 1) if (o.run.integrator is Integrator.md) else True)
    def nose_hoover_chain_length(self): pass

    @decimal
    def time(self): pass

    @decimal
    def temperature(self): pass


