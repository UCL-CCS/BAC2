from bac.utils.decorators import decimal


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

        self.couple_type_initial_lambda = kwargs.get('couple_type_initial_lambda')
        self.couple_type_final_lambda = kwargs.get('couple_type_final_lambda')
        self.couple_intramolecular = kwargs.get('couple_intramolecular')

        self.output_frequency = kwargs.get('output_frequency')
        self.output_derivatives = kwargs.get('output_derivatives')
        self.output_energy = kwargs.get('output_energy')
        self.separate_output = kwargs.get('separate_output')
        self.output_histogram_size = kwargs.get('output_histogram_size')
        self.output_bin_size = kwargs.get('output_bin_size')
