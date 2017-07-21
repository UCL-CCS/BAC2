from bac.simulate.namd.temperature_controller import TemperatureController
from bac.simulate.namd.pressure_controller import PressureController
from bac.simulate.namd.non_bonded_controller import NonBondedController
from bac.simulate.namd.constraint_controller import ConstraintController
from bac.simulate.namd.boundary_condition import PeriodicBoundaryCondition

from bac.utils.decorators import (advanced_property, positive_decimal,
                                  positive_integer,  file, boolean, back_referenced)

from bac.simulate import gromacs


class Run:

    def __init__(self, **kwargs):
        self.temperature_controller = TemperatureController(**kwargs['temperature_controller'])
        self.pressure_controller = PressureController()
        self.non_bonded_controller = NonBondedController(**kwargs['non_bonded_controller'])
        self.constraints = ConstraintController()
        self.boundary_condition = PeriodicBoundaryCondition()

        # DYNAMICS

        self.number_of_steps = kwargs.get('number_of_steps')
        self.timestep = kwargs.get('timestep')
        self.first_timestep = kwargs.get('first_timestep')

        # INPUT

        self.temperature = kwargs.get('temperature')
        self.velocities = kwargs.get('velocities')
        self.binary_velocities = kwargs.get('binary_velocities')

        self.coordinates = kwargs.get('coordinates')
        self.binary_coordinates = kwargs.get('binary_coordinates')

        self.parameters = kwargs.get('parameters')

        # OUTPUT

        self.output_name = kwargs.get('output_name')
        self.binary_output = kwargs.get('binary_output')
        self.restart_name = kwargs.get('restart_name')
        self.restart_frequency = kwargs.get('restart_frequency')
        self.restart_save = kwargs.get('restart_save')
        self.binary_restart = kwargs.get('binary_restart', True)
        self.dcd_file = kwargs.get('dcd_file')
        self.dcd_frequency = kwargs.get('dcd_frequency')
        self.dcd_unit_cell = kwargs.get('dcd_unit_cell', self.boundary_condition.is_periodic)
        self.dcd_velocity_file = kwargs.get('dcd_file')
        self.dcd_velocity_frequency = kwargs.get('dcd_velocity_frequency')
        self.dcd_force_file = kwargs.get('dcd_file')
        self.dcd_force_frequency = kwargs.get('dcd_force_frequency')

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

    # Main components:

    @back_referenced
    def temperature_controller(self): pass

    @back_referenced
    def pressure_controller(self): pass

    @back_referenced
    def non_bonded_controller(self): pass

    @back_referenced
    def constraints(self): pass

    @back_referenced
    def boundary_condition(self): pass

    # Dynamics

    @positive_integer
    def number_of_steps(self): pass

    @positive_decimal(default=1.0)
    def timestep(self): pass

    @positive_integer(default=0)
    def first_timestep(self): pass

    # Input

    @positive_decimal
    def temperature(self): pass

    @file
    def velocities(self): pass

    @file
    def binary_velocities(self): pass

    @file
    def coordinates(self): pass

    @file
    def binary_coordinates(self): pass

    @file
    def parameters(self): pass

    # Output

    @boolean(default=True)
    def binary_output(self): pass

    @file(default=lambda s: s.output_name.with_suffix('.restart'))
    def restart_name(self): pass

    @positive_integer
    def restart_frequency(self): pass

    @boolean(default=False)
    def restart_save(self): pass

    @file(default=lambda s: s.output_name.with_suffix('.dcd'))
    def dcd_file(self): pass

    @file(default=lambda s: s.output_name.with_suffix('.veldcd'))
    def dcd_velocity_file(self): pass

    @file(default=lambda s: s.output_name.with_suffix('.forcedcd'))
    def dcd_force_file(self): pass

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

    @advanced_property(type=(int, float), validator=lambda x: x >= 1, default=2.0)
    def scnb(self): pass

    # Charmm

    @file
    def structure(self): pass

    @boolean(default=True)
    def parameter_type_XPLOR(self): pass

    @boolean(default=False)
    def parameter_type_CHARMM(self): pass

    # Gromacs

    @boolean(default=False)
    def gromacs(self): pass

    # Workflow, BETA!

    @property
    def input(self):

        paths = []

        if self.binary_output:
            paths += self.binary_velocities
            paths += self.binary_coordinates
        else:
            paths += self.velocities
            paths += self.coordinates

        paths += self.parameters

        return '\n'.join(paths)

    @input.setter
    def input(self, md):
        if isinstance(md, self.__class__) and md.output_name is not None:
            if md.binary_output:
                self.binary_velocities = md.output_name.with_suffix('.vel')
                self.binary_coordinates = md.output_name.with_suffix('.coor')
            else:
                self.velocities = md.output_name.with_suffix('.vel')
                self.coordinates = md.output_name.with_suffix('.coor')

            self.parameters = md.parameters

        elif isinstance(md, gromacs.Run):
            self.gromacs = True
            self.coordinates = md.output_name.with_suffix('.gro')
            self.parameters = md.output_name.with_suffix('.top')
        else:
            raise TypeError














