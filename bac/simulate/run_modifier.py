from collections import namedtuple


ensemble_descriptor = namedtuple('ensemble_descriptor', field_names=['name', 'modifier'])


def lambda_window(iterable):
    for ld in iterable:
        name = "{}_{}".format("lambda", ld)
        yield ensemble_descriptor(name, foo(ld))


def replica(arg):

    if isinstance(arg, int):
        arg = range(arg)

    for i in arg:
        def f(run):
            pass

        name = "{}_{}".format("replica", i)

        yield ensemble_descriptor(name, f)


def foo(index):
    def f(run):
        run.free_energy_controller.initial_lambda_state = index
    return f

