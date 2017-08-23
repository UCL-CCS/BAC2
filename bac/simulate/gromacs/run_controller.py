from copy import deepcopy
from enum import Enum
import uuid

from bac.utils.decorators import positive_decimal, integer, advanced_property, boolean, file, positive_integer, decimal
from bac.simulate.gromacs.temperature_controller import TemperatureController
from bac.simulate.gromacs.pressure_controller import PressureController
from bac.simulate.gromacs.non_bonded_controller import NonBondedController
from bac.simulate.gromacs.constraint_controller import ConstraintController

from bac.simulate.simulation import Simulation


class Engine(Enum):
    namd = 'namd'
    gromacs = 'gromacs'


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

        # Minimization

        self.minimization_tolerance = kwargs.get('minimization_tolerance')
        self.minimization_step_size = kwargs.get('minimization_step_size')
        self.minimization_steepest_descent_frequency = kwargs.get('minimization_steepest_descent_frequency')
        self.minimization_correction_steps = kwargs.get('minimization_correction_steps')

        # Output

        self.coordinate_output_frequency = kwargs.get('coordinate_output_frequency')
        self.velocity_output_frequency = kwargs.get('velocity_output_frequency')
        self.force_output_frequency = kwargs.get('force_output_frequency')
        self.energy_log_output_frequency = kwargs.get('energy_log_output_frequency')
        self.recalculate_energies_frequency = kwargs.get('recalculate_energies_frequency')
        self.energy_output_frequency = kwargs.get('energy_output_frequency')
        self.compressed_coordinate_output_frequency = kwargs.get('compressed_coordinate_output_frequency')
        self.compression_precision = kwargs.get('compression_precision')
        self.compressed_groups = kwargs.get('compressed_groups')

    @advanced_property(type=Integrator, default=Integrator.md, available='5.0')
    def integrator(self): pass

    @positive_decimal(default=0, available='5.0')
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
    def preprocess_executable(self):
        """These are tasks that run relatively fast.

        Returns
        -------

        String of the executable command that you would run in bash.

        Examples
        --------

        'gmx grompp -f ...'

        """

        return 'gmx grompp -f {} -c {} -p {} -o {} -po {} {}'.\
            format(self.name.with_suffix('.mdp'),
                   self.coordinates,
                   self.topology,
                   self.output_name.with_suffix('.tpr'),
                   self.name.with_name(self.name.name + '_generated').with_suffix('.mdp'),
                   '-t' + str(self.velocities) if self.velocities else "")

    @property
    def executable(self):
        return 'gmx mdrun -nt 1 -deffnm {}'.format(self.output_name)

    def restructure_paths_with_prefix(self, prefix):
        if not self.name.is_absolute():
            self.name = prefix/self.name
        if not self.coordinates.is_absolute():
            self.coordinates = prefix/self.coordinates
        if not self.topology.is_absolute():
            self.topology = prefix/self.topology

        if self.velocities and not self.velocities.is_absolute():
            self.velocities = prefix/self.velocities

    @integer(default=0)
    def coordinate_output_frequency(self): pass

    @integer(default=0)
    def velocity_output_frequency(self): pass

    @integer(default=0)
    def force_output_frequency(self): pass

    @integer(default=1000)
    def energy_log_output_frequency(self): pass

    @integer(default=100)
    def recalculate_energies_frequency(self): pass

    @integer(default=1000, validator=lambda x, s: x % s.recalculate_energies_frequency == 0)
    def energy_output_frequency(self): pass

    @integer(default=0)
    def compressed_coordinate_output_frequency(self): pass

    @decimal(default=1000)
    def compression_precision(self): pass

    @advanced_property(type=list, default=[])
    def compressed_groups(self): pass

    @property
    def configuration_file_suffix(self):
        return ".mdp"

    @property
    def engine_type(self):
        return Engine.gromacs







