from copy import deepcopy
from pathlib import Path

import numpy as np

from bac.simulate.namd.temperature_controller import TemperatureController
from bac.simulate.namd.pressure_controller import PressureController
from bac.simulate.namd.non_bonded_controller import NonBondedController
from bac.simulate.namd.constraint_controller import ConstraintController
from bac.simulate.namd.free_energy_controller import FreeEnergyController

from bac.utils.decorators import (advanced_property, positive_decimal,
                                  positive_integer, pathlike, boolean,
                                  non_negative_integer)

from bac.simulate.basesimulation import BaseSimulation, Engine
from bac.simulate.coding import Encodable
from bac.simulate.namd.integrator import Integrator, VerletIntegrator


class Simulation(BaseSimulation, Encodable):

    def __init__(self, **kwargs):
        self.temperature_controller = None
        self.pressure_controller = None
        self.non_bonded_controller = None
        self.constraints = ConstraintController()

        # DYNAMICS

        self.integrator: Integrator = kwargs.get('integrator')
        self.number_of_steps = kwargs.get('number_of_steps')

        self.minimization = kwargs.get('minimization')

        # INPUT

        self.temperature = kwargs.get('temperature')
        self.velocities = kwargs.get('velocities')
        self.binary_velocities = kwargs.get('binary_velocities')

        self.coordinates = kwargs.get('coordinates')
        self.binary_coordinates = kwargs.get('binary_coordinates')

        self.parameters = kwargs.get('parameters')

        # OUTPUT

        self.name = kwargs.get('name')
        self.output_name = kwargs.get('output_name')
        self.binary_output = kwargs.get('binary_output')
        self.restart_name = kwargs.get('restart_name')
        self.restart_frequency = kwargs.get('restart_frequency')
        self.restart_save = kwargs.get('restart_save')
        self.binary_restart = kwargs.get('binary_restart', True)
        self.dcd_file = kwargs.get('dcd_file')
        self.dcd_frequency = kwargs.get('dcd_frequency')
        self.dcd_unit_cell = kwargs.get('dcd_unit_cell')
        self.dcd_velocity_file = kwargs.get('dcd_file')
        self.dcd_velocity_frequency = kwargs.get('dcd_velocity_frequency')
        self.dcd_force_file = kwargs.get('dcd_file')
        self.dcd_force_frequency = kwargs.get('dcd_force_frequency')

        self.output_energies = kwargs.get('output_energies')
        self.output_pressure = kwargs.get('output_energies')

        # OTHER SETTINGS

        self.center_of_mass_motion = kwargs.get('center_of_mass_motion')
        self.random_seed = kwargs.get('random_seed')
        self.zero_momentum = kwargs.get('zero_momentum')

        # Forcefield specifics

        # AMBER
        self.amber = kwargs.get('amber')
        self.read_exclusions = kwargs.get('read_exclusions')
        self.scnb = kwargs.get('scnb')

        # CHARMM

        self.structure = kwargs.get('structure')
        self.parameter_type_XPLOR = kwargs.get('parameter_type_XPLOR')
        self.parameter_type_CHARMM = kwargs.get('parameter_type_CHARMM')

        # GROMACS

        self.gromacs = kwargs.get('gromacs')

    # Dynamics

    @advanced_property(type=Integrator, default=VerletIntegrator())
    def integrator(self): pass

    @positive_integer
    def number_of_steps(self): pass

    @boolean(default=False)
    def minimization(self): pass

    # Input

    @positive_decimal
    def temperature(self): pass

    @pathlike
    def velocities(self): pass

    @velocities.post_set_processing
    def velocities(self, value):
        if value is not None:
            self.temperature = None

    @pathlike
    def binary_velocities(self): pass

    @pathlike
    def coordinates(self): pass

    @pathlike
    def binary_coordinates(self): pass

    @pathlike
    def parameters(self): pass

    # Output

    @pathlike
    def name(self): pass

    @pathlike
    def output_name(self): pass

    @boolean(default=True)
    def binary_output(self): pass

    @pathlike(default=lambda s: s.output_name.with_suffix('.restart'))
    def restart_name(self): pass

    @positive_integer
    def restart_frequency(self): pass

    @boolean(default=False)
    def restart_save(self): pass

    @pathlike(default=lambda s: s.output_name.with_suffix('.dcd'))
    def dcd_file(self): pass

    @pathlike(default=lambda s: s.output_name.with_suffix('.veldcd'))
    def dcd_velocity_file(self): pass

    @boolean(default=True)
    def dcd_unit_cell(self): pass

    @pathlike(default=lambda s: s.output_name.with_suffix('.forcedcd'))
    def dcd_force_file(self): pass

    @positive_integer
    def output_energies(self): pass

    @non_negative_integer
    def output_pressure(self): pass

    # Other settings

    @boolean(default=False)
    def center_of_mass_motion(self): pass

    @positive_integer
    def random_seed(self): pass

    @boolean(default=False)
    def zero_momentum(self): pass

    # Forcefield specifics

    # Amber

    @boolean(default=False)
    def amber(self): pass

    @boolean(default=True)
    def read_exclusions(self): pass

    @advanced_property(type=np.float, validator=lambda _, x: x >= 1, default=2.0)
    def scnb(self): pass

    # Charmm

    @pathlike
    def structure(self): pass

    @boolean(default=True)
    def parameter_type_XPLOR(self): pass

    @boolean(default=False)
    def parameter_type_CHARMM(self): pass

    # Gromacs

    @boolean(default=False)
    def gromacs(self): pass

    @property
    def configuration_file_suffix(self):
        return ".conf"

    @property
    def engine_type(self):
        return Engine.namd

    @property
    def executable(self):
        return f"aprun -n 1 -N 1 -d 31 /u/sciteam/jphillip/NAMD_LATEST_CRAY-XE-ugni-smp-BlueWaters/namd2" \
               f" +ppn 30 +pemap 0-29 +commap 30 {self.name.with_suffix('.conf')} >& {self.name.with_suffix('.log')}"

    @property
    def preprocess_executable(self):
        return ''

    def add_input_dependency(self, other_simulation):
        super(Simulation, self).add_input_dependency(other_simulation)

        self.non_bonded_controller.extended_system = other_simulation.output_name.with_suffix('.xsc')
        self.coordinates = other_simulation.output_name.with_suffix('.coor')

        if self.constraints.harmonic_constraint is not None:
            self.constraints.harmonic_constraint.reference_position_file = other_simulation.output_name.with_suffix('.coor')

        if not other_simulation.minimization:
            self.velocities = other_simulation.output_name.with_suffix('.vel')

    def restructure_paths_with_prefix(self, prefix):

        if not self.name.exists():
            self.name = prefix/self.name

        prefix = Path('/'.join(['..']*len(prefix.parts)))

        if self.coordinates and self.coordinates.exists():
            self.coordinates = prefix/self.coordinates

        if self.parameters and self.parameters.exists():
            self.parameters = prefix/self.parameters

        if self.non_bonded_controller.extended_system and self.non_bonded_controller.extended_system.exists():
            self.non_bonded_controller.extended_system = prefix/self.non_bonded_controller.extended_system

        if self.non_bonded_controller.xst_file and self.non_bonded_controller.xst_file.exists():
            self.non_bonded_controller.xst_file = prefix/self.non_bonded_controller.xst_file

        if self.constraints.harmonic_constraint and self.constraints.harmonic_constraint.reference_position_file and self.constraints.harmonic_constraint.reference_position_file.exists():
            self.constraints.harmonic_constraint.reference_position_file = prefix/self.constraints.harmonic_constraint.reference_position_file

        if self.constraints.harmonic_constraint and self.constraints.harmonic_constraint.force_constant_file and self.constraints.harmonic_constraint.force_constant_file.exists():
            self.constraints.harmonic_constraint.force_constant_file = prefix/self.constraints.harmonic_constraint.force_constant_file

        if self.free_energy_controller and self.free_energy_controller.file:
            self.free_energy_controller.file = prefix/self.free_energy_controller.file


    def __next__(self):
        next_run = deepcopy(self)
        next_run.name, next_run.output_name = None, None
        next_run.add_input_dependency(self)
        return next_run

    def encode(self, path=None, suffix=None):
        path = path if (path is None or path.suffix == self.configuration_file_suffix) else path / self.name.with_suffix(
            self.configuration_file_suffix)
        suf = (suffix or '') + f"{'minimize' if self.minimization else 'run'} {self.number_of_steps}\n"
        return super(Simulation, self).encode(path=path, suffix=suf)


















