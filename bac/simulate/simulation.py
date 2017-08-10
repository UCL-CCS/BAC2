from abc import ABCMeta, abstractmethod


class Simulation(metaclass=ABCMeta):

    @property
    @abstractmethod
    def executable(self):
        return NotImplemented

    @property
    @abstractmethod
    def preprocess_executable(self):
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
        return self.dependency.is_finished and self.state < State.running

    def restructure_paths_with_prefix(self, prefix):
        return NotImplemented