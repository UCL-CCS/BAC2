from collections import namedtuple


ensemble_descriptor = namedtuple('ensemble_descriptor', field_names=['name', 'modifier'])


def lambda_window(arg):
    if isinstance(arg, int):
        arg = range(arg)
    for ld in arg:
        name = "{}_{}".format("lambda", ld)
        yield ensemble_descriptor(name, _lambda_window(ld))


def _lambda_window(index):
    # This wrapper is needed around the actual function that is yielded by `lambda_window`
    # because otherwise `index` is not captured but referenced from the outside scope,
    # meaning it will be the last value by the time this function is actually called.
    def f(run):
        run.free_energy_controller.initial_lambda_state = index
    return f


def replica(arg):
    """Replica modifier

    :param arg: iterable or int. If int then converted to range(int)
    :return: ensemble_descriptor
    """
    if isinstance(arg, int):
        arg = range(arg)
    for i in arg:
        def f(run): pass
        name = "{}_{}".format("replica", i)
        yield ensemble_descriptor(name, f)


