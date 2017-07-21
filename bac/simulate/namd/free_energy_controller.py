from enum import Enum

from bac.utils.decorators import advanced_property, decimal, positive_integer, file, column, positive_decimal, boolean


class FreeEnergyCalculationType(Enum):
    fep = 'fep'
    ti = 'ti'


class FreeEnergyController:

    def __init__(self, **kwargs):
        self.type = kwargs.get('type')
        self.lambda_start = kwargs.get('lambda_start')
        self.lambda_end = kwargs.get('lambda_end')
        self.equilibration_steps = kwargs.get('equilibration_steps')
        self.file = kwargs.get('file')
        self.column = kwargs.get('column')
        self.output_frequency = kwargs.get('output_frequency')
        self.output_name = kwargs.get('output_name')
        self.van_der_waals_shift_coefficient = kwargs.get('van_der_waals_shift_coefficient')
        self.electronic_interaction_start_lambda = kwargs.get('electronic_interaction_start_lambda')
        self.van_der_waals_end_lambda = kwargs.get('van_der_waals_end_lambda')
        self.bonded_interaction_end_lambda = kwargs.get('bonded_interaction_end_lambda')
        self.bond_decoupling = kwargs.get('bond_decoupling')
        self.decouple = kwargs.get('decouple')

    @advanced_property(type=FreeEnergyCalculationType, default=FreeEnergyCalculationType.ti)
    def type(self): pass

    @decimal(validator=lambda x: 0 <= x <= 1)
    def lambda_start(self): pass

    @decimal(validator=lambda x: 0 <= x <= 1)
    def lambda_end(self): pass

    @positive_integer(default=0, validator=lambda x, s: s.run.number_of_steps >= x)
    def equilibration_steps(self): pass

    @file(default=lambda s: s.run.coordinates)
    def file(self): pass

    @column
    def column(self): pass

    @positive_integer(default=5)
    def output_frequency(self): pass

    @advanced_property(type=str, default='outfilename')
    def output_name(self): pass

    @positive_decimal(default=5)
    def van_der_waals_shift_coefficient(self): pass

    @positive_decimal(default=0.5)
    def electronic_interaction_start_lambda(self): pass

    @positive_decimal(default=1)
    def van_der_waals_end_lambda(self): pass

    @positive_decimal(default=0.0)
    def bonded_interaction_end_lambda(self): pass

    @boolean(default=False)
    def bond_decoupling(self): pass

    @boolean(default=False)
    def decouple(self): pass
    





