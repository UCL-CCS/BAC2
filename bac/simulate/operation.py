from copy import deepcopy
from enum import Enum
import uuid
from bac.simulate import namd, gromacs


class OperationState(Enum):
    prep = 0
    ready = 1
    executing = 2
    finished = 3
    failed = 4


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

    @property
    def execute_path(self):
        args_list = []

        if isinstance(self.run, gromacs.Run):
            args1 = ['gmx', 'grompp',
                     '-f', '{}.mpd'.format(uuid.uuid4()),
                     '-c', str(self.run.coordinates),
                     '-p', str(self.run.topology),
                     '-o', str(self.run.output_name.with_suffix('.tpr'))]
            args1 += [ '-t', str(self.run.velocities)] if self.run.velocities is not None else []

            args2 = ['gmx', 'mdrun', '-nt', 1, '-deffnm', str(self.run.output_name)]

            args_list = [args1, args2]

        elif isinstance(self.run, namd.Run):
            args = ['namd2', 'input_file_name.conf', '>', self.run.output_name.with_suffix('.log')]
            args_list = [args]

        return args_list