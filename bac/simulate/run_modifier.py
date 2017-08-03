class ReplicaModifier:

    def __init__(self, count, name='replica'):
        self.index = count
        self.name = name

    def __iter__(self):
        return self

    def __next__(self):
        if self.index == 0:
            raise StopIteration
        self.index = self.index - 1
        return self.do_nothing, '{}_{}'.format(self.name, self.index)

    @staticmethod
    def do_nothing(run):
        pass


class LambdaModifier:
    def __init__(self, lambda_window_count, name='lambda'):
        self.index = lambda_window_count
        self.name = name

    def __iter__(self):
        return self

    def __next__(self):
        if self.index == 0:
            raise StopIteration
        self.index = self.index - 1
        return self.lambda_window(self.index), '{}_{}'.format(self.name, self.index)

    @staticmethod
    def lambda_window(index):
        def func(run):
            run.free_energy_controller.initial_lambda_state = index
        return func
