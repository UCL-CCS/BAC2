class OperationQueue:

    def __init__(self, resource):
        self.resource = resource
        self.operations = []

    def add_operation(self, operation):
        self.operations.append(operation)

    def execute(self):
        print('Executing on {}'.format(self.resource))



