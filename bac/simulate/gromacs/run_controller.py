from enum import Enum

from bac.utils.decorators import positive_decimal, integer, advanced_property, boolean


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


class CenterOfMassMotion(Enum):
    linear = 'Linear'
    angular = 'Angular'
    none = 'None'


class Run:

    def __init__(self, **kwargs):
        self.integrator = kwargs.get('integrator')
        self.initial_time = kwargs.get('initial_time')
        self.delta_time = kwargs.get('delta_time')
        self.number_of_steps = kwargs.get('number_of_steps')
        self.initial_step = kwargs.get('initial_step')
        self.center_of_mass_motion = kwargs.get('center_of_mass_motion')
        self.com_motion_removal = kwargs.get('com_motion_removal')
        self.comm_groups = kwargs.get('comm_groups')

        self.generate_velocities = kwargs.get('generate_velocities')
        self.generate_temperature = kwargs.get('generate_temperature')
        self.generate_seed = kwargs.get('generate_seed')

        self.minimization_tolerance = None
        self.minimization_step_size = None

    @advanced_property(type=Integrator, default=Integrator.md)
    def integrator(self): pass

    @positive_decimal(default=0)
    def initial_time(self): pass

    @positive_decimal(default=0.001)
    def delta_time(self): pass

    @integer(default=0)
    def number_of_steps(self): pass

    @integer(default=0)
    def initial_step(self): pass

    @advanced_property(type=CenterOfMassMotion, default=CenterOfMassMotion.linear)
    def center_of_mass_motion(self): pass

    @integer(default=100)
    def com_motion_removal(self): pass

    @boolean(default=False)
    def generate_veolcities(self): pass

    @positive_decimal(default=300)
    def generate_temperature(self): pass

    @integer(default=-1)
    def generate_seed(self): pass






