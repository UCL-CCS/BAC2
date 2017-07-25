from .operation import Operation


class OperationQueue:

    def __init__(self, resource):
        self.resource = resource
        self.operations = []

    def add_operation(self, operation):
        self.operations.append(operation)

    def add_run(self, run):
        op = Operation(run=run)
        self.add_operation(op)

    def execute(self):
        print('Executing on {}'.format(self.resource))



