import numpy as np

from _collections_abc import Iterable, Iterator
from typing import Callable
from .basesimulation import BaseSimulation


class EnsemblePackage:
    def __init__(self, path_name: str=None, fn: Callable[[BaseSimulation], None]=None, iterator_state=None):
        self.path_name = path_name
        self.fn = fn or self._do_nothing
        self.iterator_state = iterator_state

    def _do_nothing(self, *args):
        pass


class BaseEnsembleIterator(Iterator):
    def __init__(self, name: str, underlying_iterable: Iterable):
        self.name = name
        self.underlying_iterable = underlying_iterable

    @property
    def generic_bash_name(self):
        return f'${{{self.name.upper()}}}'

    @property
    def generic_path_name(self):
        return f'{self.name}_{self.generic_bash_name}'

    def __iter__(self):
        self._iterator = iter(self.underlying_iterable)
        return self

    def __next__(self) -> EnsemblePackage:
        ens_pack = EnsemblePackage()
        ens_pack.iterator_state = next(self._iterator)
        ens_pack.path_name = f'{self.name}_{ens_pack.iterator_state.__str__()}'
        return ens_pack


class Replica(BaseEnsembleIterator):
    def __init__(self, number_of: int):
        super().__init__(name='replica', underlying_iterable=range(number_of))


class LambdaWindow(BaseEnsembleIterator):
    def __init__(self, number_of_states: int=None, number_of_windows: int=None, additional=None):
        it = range(number_of_states) if number_of_states else np.linspace(0.0, 1.0, number_of_windows, endpoint=True)
        if additional:
            it = np.append(it, additional)
            it.sort()
        super().__init__(name='lambda', underlying_iterable=it)

    def __next__(self) -> EnsemblePackage:
        ens_pack = super(LambdaWindow, self).__next__()
        ens_pack.fn = _lambda_window(ens_pack.iterator_state)
        return ens_pack


def _lambda_window(n):
    def f(simulation: BaseSimulation):
        if simulation.free_energy_controller.supports_lambda_state is True:
            simulation.free_energy_controller.initial_lambda_state = n
        else:
            simulation.free_energy_controller.start_lambda = n
            simulation.free_energy_controller.end_lambda = n

    return f


