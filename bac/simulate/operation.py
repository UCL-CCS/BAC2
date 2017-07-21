from copy import deepcopy


class Operation:

    def __init__(self, run):
        self.run = deepcopy(run)
        self.dependencies = []

    def add_dependency(self, other):
        self.dependencies.append(other)