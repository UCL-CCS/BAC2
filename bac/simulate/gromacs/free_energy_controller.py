from bac.utils.decorators import decimal


class FreeEnergyController:
    def __init__(self, **kwargs):
        self.expanded = kwargs.get('expanded')
        self.initial_lambda = kwargs.get('')
        self.delta_lambdas = kwargs.get('')
        self.initial_lambda_state = kwargs.get('')
        sek