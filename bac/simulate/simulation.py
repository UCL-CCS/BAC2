from abc import ABCMeta, abstractmethod


class Simulation(metaclass=ABCMeta):

    @property
    @abstractmethod
    def executable_path(self):
        return NotImplemented

    @abstractmethod
    def add_input_dependency(self, other_simulation):
        self.dependency = other_simulation

    @property
    def dependency(self):
        return self._dependency

    @dependency.setter
    def dependency(self, v):
        self._dependency = v

    def is_ready(self):
        return self.dependency.is_finished