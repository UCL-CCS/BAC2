from enum import Enum

from bac.utils.decorators import decimal, advanced_property, boolean, integer


class CouplingType(Enum):
    vdw_and_coulomb = 'vdw-q'
    vdw = 'vdw'
    coulomb = 'q'
    none = 'none'


class OutputEnergyType(Enum):
    no = 'no'
    potential = 'potential'
    total = 'total'


class FreeEnergyController:

    def __init__(self, **kwargs):
        self.expanded = kwargs.get('expanded')
        self.initial_lambda = kwargs.get('initial_lambda')
        self.delta_lambda = kwargs.get('delta_lambda')
        self.initial_lambda_state = kwargs.get('initial_lambda_state')
        self.fep_lambdas = kwargs.get('fep_lambdas')
        self.coulomb_lambdas = kwargs.get('coulomb_lambdas')
        self.van_der_waals_lambdas = kwargs.get('van_der_waals_lambdas')
        self.bonded_lambdas = kwargs.get('bonded_lambdas')
        self.restraint_lambdas = kwargs.get('restraint_lambdas')
        self.mass_lambdas = kwargs.get('mass_lambdas')
        self.temperature_lambdas = kwargs.get('temperature_lambdas')
        self.calculate_lambda_neighbors = kwargs.get('calculate_lambda_neighbors')

        self.soft_core_alpha = kwargs.get('soft_core_alpha')
        self.soft_core_radial_power = kwargs.get('soft_core_radial_power')
        self.soft_core_coulomb = kwargs.get('soft_core_coulomb')
        self.soft_core_power = kwargs.get('soft_core_power')
        self.soft_core_sigma = kwargs.get('soft_core_sigma')

        self.couple_molecule_groups = kwargs.get('couple_molecule_groups')
        self.couple_type_initial_lambda = kwargs.get('couple_type_initial_lambda')
        self.couple_type_final_lambda = kwargs.get('couple_type_final_lambda')
        self.couple_intramolecular = kwargs.get('couple_intramolecular')

        self.output_frequency = kwargs.get('output_frequency')
        self.output_derivatives = kwargs.get('output_derivatives')
        self.output_energy = kwargs.get('output_energy')
        self.separate_output = kwargs.get('separate_output')
        self.output_histogram_size = kwargs.get('output_histogram_size')
        self.output_bin_size = kwargs.get('output_bin_size')

    @boolean(default=False)
    def expanded(self): pass

    @decimal(default=-1)
    def initial_lambda(self): pass

    @decimal(default=0)
    def delta_lambda(self): pass

    @integer(default=-1)
    def initial_lambda_state(self): pass

    @advanced_property(type=list, default=[])
    def fep_lambdas(self): pass

    @advanced_property(type=list, default=[])
    def coulomb_lambdas(self): pass

    @advanced_property(type=list, default=[])
    def van_der_waals_lambdas(self): pass

    @advanced_property(type=list, default=[])
    def bonded_lambdas(self): pass

    @advanced_property(type=list, default=[])
    def restraint_lambdas(self): pass

    @advanced_property(type=list, default=[])
    def mass_lambdas(self): pass

    @advanced_property(type=list, default=[])
    def temperature_lambdas(self): pass

    @integer(default=1)
    def calculate_lambda_neighbors(self): pass

    @decimal(default=0)
    def soft_core_alpha(self): pass

    @integer(default=6, validator=lambda x, _: x in [6, 48])
    def soft_core_radial_power(self): pass

    @boolean(default=False)
    def soft_core_coulomb(self): pass

    @integer(default=0, validator=lambda x, _: x in [1, 2])
    def soft_core_power(self): pass

    @decimal(default=0.3)
    def soft_core_sigma(self): pass

    @advanced_property(type=list, default=[])
    def couple_molecule_groups(self): pass

    @advanced_property(type=CouplingType, default=CouplingType.vdw_and_coulomb)
    def couple_type_initial_lambda(self): pass

    @advanced_property(type=CouplingType, default=CouplingType.vdw_and_coulomb)
    def couple_type_final_lambda(self): pass

    @boolean(default=False)
    def couple_intramolecular(self): pass

    @integer(default=100)
    def output_frequency(self): pass

    @boolean(default=True)
    def output_derivatives(self): pass

    @advanced_property(type=OutputEnergyType, default=OutputEnergyType.no)
    def output_energy(self): pass

    @boolean(default=True)
    def separate_output(self): pass

    @decimal(default=0)
    def output_histogram_size(self): pass

    @decimal(default=0.1)
    def output_bin_size(self): pass