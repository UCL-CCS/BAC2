from copy import deepcopy
from enum import Enum

from bac.utils.decorators import positive_decimal, integer, advanced_property, boolean, file, back_referenced, positive_integer
from bac.simulate.gromacs.temperature_controller import TemperatureController
from bac.simulate.gromacs.pressure_controller import PressureController
from bac.simulate.gromacs.non_bonded_controller import NonBondedController
from bac.simulate.gromacs.constraint_controller import ConstraintController


from bac.simulate import namd


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

        self.temperature_controller = TemperatureController()
        self.pressure_controller = PressureController()
        self.non_bonded_controller = NonBondedController()
        self.constraints = ConstraintController()
        self.free_energy_controller = None

        self.integrator = kwargs.get('integrator')
        self.initial_time = kwargs.get('initial_time')
        self.delta_time = kwargs.get('delta_time')
        self.number_of_steps = kwargs.get('number_of_steps')
        self.initial_step = kwargs.get('initial_step')
        self.center_of_mass_motion = kwargs.get('center_of_mass_motion')
        self.com_motion_removal_frequency = kwargs.get('com_motion_removal_frequency')
        self.comm_groups = kwargs.get('comm_groups')

        self.generate_velocities = kwargs.get('generate_velocities')
        self.generate_temperature = kwargs.get('generate_temperature')
        self.generate_seed = kwargs.get('generate_seed')

        self.coordinates = kwargs.get('coordinates')
        self.topology = kwargs.get('topology')
        self.velocities = kwargs.get('velocities')

        self.output_name = kwargs.get('output_name')

        self.minimization_tolerance = None
        self.minimization_step_size = None
        self.minimization_steepest_descent_frequency = None
        self.minimization_correction_steps = None

    # Main components:

    @back_referenced
    def temperature_controller(self): pass

    @back_referenced
    def pressure_controller(self): pass

    @back_referenced
    def non_bonded_controller(self): pass

    @back_referenced
    def constraints(self): pass

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
    def com_motion_removal_frequency(self): pass

    @advanced_property(type=list, default=['system'])
    def comm_groups(self): pass

    @boolean(default=False)
    def generate_velocities(self): pass

    @positive_decimal(default=300)
    def generate_temperature(self): pass

    @integer(default=-1)
    def generate_seed(self): pass

    @positive_decimal(default=10.0)
    def minimization_tolerance(self): pass

    @positive_decimal(default=0.01)
    def minimization_step_size(self): pass

    @positive_integer(default=1000)
    def minimization_steepest_descent_frequency(self): pass

    @positive_integer(default=10)
    def minimization_correction_steps(self): pass

    @file
    def coordinates(self): pass

    @file
    def topology(self): pass

    @file
    def velocities(self): pass

    @file
    def output_name(self): pass

    @property
    def input(self):

        paths = []
        paths += self.velocities
        paths += self.coordinates
        paths += self.topology

        return '\n'.join(paths)

    @input.setter
    def input(self, md):
        if isinstance(md, self.__class__):
            if md.output_name is None: raise ValueError('Output not defined')
            self.coordinates = md.output_name.with_suffix('.gro')
            self.topology = md.topology
            self.velocities = md.output_name.with_suffix('.cpt') if md.constraints.continuation is True else None
        elif isinstance(md, namd.Run):
            raise NotImplementedError
        else:
            raise TypeError

    def __iter__(self):
        return self

    def __next__(self):
        next_run = deepcopy(self)
        next_run.output_name = None
        next_run.coordinates = self.output_name.with_suffix('.gro')
        next_run.topology = self.topology
        next_run.velocities = self.output_name.with_suffix('.cpt') if self.constraints.continuation is True else None
        return next_run






