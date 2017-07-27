from copy import deepcopy
from enum import Enum


class OperationState(Enum):
    prep = 0
    ready = 1
    executing = 2
    finished = 3
    cancelled = 4

class Operation:

    def __init__(self, run):
        self.run = deepcopy(run)
        self.dependencies: [Operation] = []
        self.state = OperationState.prep

    def add_dependency(self, other):
        self.dependencies.append(other)

    @property
    def state(self):
        return self._state

    @state.setter
    def state(self, new_state):
        self._state = new_state

    @property
    def is_ready(self):
        for dependency in self.dependencies:
            if not (dependency.state is OperationState.finished):
                return False
        return True
