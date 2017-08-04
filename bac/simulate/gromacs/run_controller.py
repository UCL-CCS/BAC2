from copy import deepcopy
from enum import Enum
import uuid

from bac.utils.decorators import positive_decimal, integer, advanced_property, boolean, file, back_referenced, positive_integer
from bac.simulate.gromacs.temperature_controller import TemperatureController
from bac.simulate.gromacs.pressure_controller import PressureController
from bac.simulate.gromacs.non_bonded_controller import NonBondedController
from bac.simulate.gromacs.constraint_controller import ConstraintController

from bac.simulate.simulation import Simulation


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


class Run(Simulation):

    def __init__(self, **kwargs):

        # Required controllers

        self.temperature_controller = TemperatureController()
        self.pressure_controller = PressureController()
        self.non_bonded_controller = NonBondedController()
        self.constraints = ConstraintController()

        # Optional controllers

        self.free_energy_controller = None

        # Main attributes

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
        self.name = kwargs.get('name')

        self.minimization_tolerance = kwargs.get('minimization_tolerance')
        self.minimization_step_size = kwargs.get('minimization_step_size')
        self.minimization_steepest_descent_frequency = kwargs.get('minimization_steepest_descent_frequency')
        self.minimization_correction_steps = kwargs.get('minimization_correction_steps')

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

    @file(validator=lambda x, _: x.suffix == '.gro', warn=True)
    def coordinates(self): pass

    @file(validator=lambda x, _: x.exists())
    def topology(self): pass

    @file
    def velocities(self): pass

    @file(default=lambda s: s.name)
    def output_name(self): pass

    @file(default=lambda s: s.integrator.name)
    def name(self): pass

    def add_input_dependency(self, other_simulation):
        super(Run, self).add_input_dependency(other_simulation)

        self.coordinates = other_simulation.output_name.with_suffix('.gro')
        self.topology = other_simulation.topology
        self.velocities = other_simulation.output_name.with_suffix('.cpt') if other_simulation.constraints.continuation is True else None

    def __next__(self):
        next_run = deepcopy(self)
        next_run.name, next_run.output_name = None, None
        next_run.add_input_dependency(self)
        return next_run

    @property
    def executable_path(self):
        args1 = ['gmx', 'grompp',
                 '-f', str(self.name.with_suffix('.mdp')),
                 '-c', str(self.coordinates),
                 '-p', str(self.topology),
                 '-o', str(self.output_name.with_suffix('.tpr'))]
        args1 += ['-t', str(self.velocities)] if self.velocities is not None else []

        args2 = ['gmx', 'mdrun', '-nt', 1, '-deffnm', str(self.output_name)]

        return [args1, args2]







