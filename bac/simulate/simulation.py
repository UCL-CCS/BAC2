from abc import ABCMeta, abstractmethod
from enum import Enum


class SimulationState(Enum):
    pending = 0
    ready = 1
    executing = 2
    finished = 3
    cancelled = 4

    def __eq__(self, other):
        return self.value == other.value

    def __lt__(self, other):
        return self.value < other.value

    def __le__(self, other):
        return self.value <= other.value


class Simulation(Versioned, metaclass=ABCMeta):

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
        return self.dependency.is_finished and self.state < SimulationState.executing

    def restructure_paths_with_prefix(self, prefix):
        return NotImplemented
