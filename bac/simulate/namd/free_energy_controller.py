import numpy as np

from enum import Enum

from bac.utils.decorators import advanced_property,  positive_integer, pathlike, pdbcolumn, positive_decimal, boolean

from bac.simulate.coding import Encodable


class FreeEnergyCalculationType(Enum):
    fep = 'fep'
    ti = 'ti'


class FreeEnergyController(Encodable):

    def __init__(self, **kwargs):
        self.type = kwargs.get('type')
        self.start_lambda = kwargs.get('start_lambda')
        self.end_lambda = kwargs.get('end_lambda')
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

    @advanced_property(type=np.float, validator=lambda _, x: 0 <= x <= 1)
    def start_lambda(self): pass

    @advanced_property(type=np.float, validator=lambda _, x: 0 <= x <= 1)
    def end_lambda(self): pass

    @property
    def supports_lambda_state(self):
        return False

    @positive_integer(default=0, validator=lambda self, x: self.run.number_of_steps >= x)
    def equilibration_steps(self): pass

    @pathlike(default=lambda self: self.simulation.coordinates)
    def file(self): pass

    @pdbcolumn
    def column(self): pass

    @positive_integer(default=5)
    def output_frequency(self): pass

    @pathlike(default=lambda self:self.simulation.output_name.with_suffix('.alch'))
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
    





