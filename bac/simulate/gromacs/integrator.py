from supproperty import positive_integer, positive_decimal
from bac.simulate.coding import Encodable


class Integrator:
    def __init__(self, timestep=None, first_step=None, initial_time=None):
        self.timestep = timestep
        self.first_step = first_step
        self.initial_time = initial_time

    @positive_decimal(default=0.001)
    def timestep(self): pass

    @positive_integer(default=0)
    def first_step(self): pass

    @positive_decimal(default=0, available='5.0')
    def initial_time(self): pass


class VerletIntegrator(Integrator, Encodable):
    def __init__(self, *, timestep=None, first_step=None, initial_time=None):
        super(VerletIntegrator, self).__init__(timestep, first_step, initial_time)


class VelocityVerletIntegrator(Integrator, Encodable):
    def __init__(self, *, timestep=None, first_step=None, initial_time=None):
        super(VelocityVerletIntegrator, self).__init__(timestep, first_step, initial_time)


class AverageVelocityVerletIntegrator(Integrator, Encodable):
    def __init__(self, *, timestep=None, first_step=None, initial_time=None):
        super(AverageVelocityVerletIntegrator, self).__init__(timestep, first_step, initial_time)


class LangevinIntegrator(Integrator, Encodable):
    def __init__(self, *, timestep=None, first_step=None, initial_time=None, temperature, damping=None):
        super(LangevinIntegrator, self).__init__(timestep, first_step, initial_time)
        self.temperature = temperature
        self.damping = damping

    @positive_decimal
    def temperature(self): pass

    @positive_decimal
    def damping(self): pass


# class Integrator(Enum):
#     md = 'md'
#     md_vv = 'md-vv'
#     md_vv_avek = 'md-vv-avek'
#     sd = 'sd'
#     bd = 'bd'
#     steep = 'steep'
#     cg = 'cg'
#     l_bfgs = 'l-bfgs'
#     nm = 'nm'
#     tpi = 'tpi'
#     tpic = 'tpic'


