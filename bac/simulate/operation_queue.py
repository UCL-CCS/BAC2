import subprocess
import uuid
import itertools
from pathlib import Path

from .operation import Operation


class OperationQueue:

    def __init__(self, resource, name):
        self.resource = resource
        self.name = name or str(uuid.uuid4())
        self.path = Path(self.name)
        self.operations = []
        self.modifiers = None

    def add_operation(self, operation):
        self.operations.append(operation)

    def add_run(self, run):
        op = Operation(run=run)
        self.add_operation(op)

    def execute(self):

        while len(self.operations) > 0:
            op = next(op for op in self.operations if op.is_ready)

            subprocess.run(op.execute_path)

        print('Executing on {}'.format(self.resource))

    def create_file_structure(self):
        tree = [['{}_{}'.format(mod.name, s) for s in mod] for mod in self.modifiers]

        for branch in itertools.product(*tree):
            self.path.joinpath(*branch).mkdir(parents=True)






