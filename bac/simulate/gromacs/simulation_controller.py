from copy import deepcopy
from enum import Enum
from pathlib import Path

from bac.utils.decorators import positive_decimal, integer, advanced_property, boolean, pathlike, positive_integer, decimal
from bac.simulate.gromacs.temperature_controller import TemperatureController
from bac.simulate.gromacs.pressure_controller import PressureController
from bac.simulate.gromacs.non_bonded_controller import NonBondedController
from bac.simulate.gromacs.constraint_controller import ConstraintController

from bac.simulate.basesimulation import BaseSimulation, Engine
from bac.simulate.gromacs.integrator import Integrator
from bac.simulate.coding import Encodable


class CenterOfMassMotion(Enum):
    linear = 'Linear'
    angular = 'Angular'
    none = 'None'

    @classmethod
    def _missing_(cls, value):
        if value is None:
            return cls.none
        else:
            super()._missing_(value)


class Simulation(BaseSimulation, Encodable):

    def __init__(self, **kwargs):

        # Required controllers

        self.temperature_controller = TemperatureController()
        self.pressure_controller = PressureController()
        self.non_bonded_controller = NonBondedController()
        self.constraints = ConstraintController()

        # Optional controllers

        # Main attributes

        self.integrator: Integrator = kwargs.get('integrator')
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

        self.coordinates: Path = kwargs.get('coordinates')
        self.topology: Path = kwargs.get('topology')
        self.velocities: Path = kwargs.get('velocities')

        self.output_name: Path = kwargs.get('output_name')
        self.name: Path = kwargs.get('name')

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

    @pathlike(validator=lambda _, x: x.suffix == '.gro', warn=True)
    def coordinates(self): pass

    @pathlike(validator=lambda _, x: x.exists())
    def topology(self): pass

    @pathlike
    def velocities(self): pass

    @pathlike(default=lambda s: s.name)
    def output_name(self): pass

    @pathlike(default=lambda s: s.integrator.name)
    def name(self): pass

    def add_input_dependency(self, other_simulation):
        super(Simulation, self).add_input_dependency(other_simulation)

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

        return f"aprun -n 1 gmx_mpi grompp -f {self.name.with_suffix('.mdp')} -c {self.coordinates} " \
               f"-p {self.topology} -o {self.output_name.with_suffix('.tpr')} " \
               f"-po {self.name.with_name(self.name.name + '_generated').with_suffix('.mdp')} " \
               f"{'-t ' + str(self.velocities) if self.velocities else ''} " \
               f"&> {self.name.with_name(self.name.name + '_grompp').with_suffix('.out')}"

    @property
    def executable(self):
        return f"aprun -n {1 if self.integrator is Integrator.steep else 32} gmx_mpi mdrun -deffnm {self.output_name} " \
               f"&> {self.output_name.with_suffix('.out')}"

    def restructure_paths_with_prefix(self, prefix):
        if not self.name.exists():
            self.name = prefix/self.name
        if not self.coordinates.exists():
            self.coordinates = prefix/self.coordinates
        if self.velocities and not self.velocities.exists():
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

    @integer(default=1000, validator=lambda self, x: x % self.recalculate_energies_frequency == 0)
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

    def encode(self, path=None, suffix=None):
        path = path if path.suffix == self.configuration_file_suffix else path / self.name.with_suffix(
            self.configuration_file_suffix)
        return super(Simulation, self).encode(path=path)

