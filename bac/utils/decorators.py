import copy
import warnings

from supproperty import supproperty
from .pdb import PDBColumn


class pdbcolumn(supproperty):
    def __init__(self, *args, **kwargs):
        kwargs['type'] = PDBColumn
        if kwargs.get('default') is None:
            kwargs['default'] = 'O'
        super(pdbcolumn, self).__init__(*args, **kwargs)


class advanced_property(property):
    """Decorator for attributes that reference their owner.

    The decorated object has a property called `simulation` that references back to the
    main `Simulation` object. This is important for example when validators have to access
    properties from other parts of the object tree.

    Examples
    --------

    @back_referenced
    def pressure_controller(self): pass

    # Then in the class definition access self.simulation:

    class PressureController:

        @decimal(validator:lambda x, s: x < s.simulation.pressure)
        def pressure(self): pass

    Notes
    -----

    If you set an @back_referenced property that already has another container
    then the object is copied before it is set as the @back_referenced. This
    means that you can only access it from the main class' property. If this
    happens a warning is issued.

    """
    def __init__(self, f):
        self.f = f
        super(back_referenced, self).__init__(fget=self._fget, fset=self._fset)

    def _fget(self, obj):
        try:
            return obj.__getattribute__(self.private_name)
        except AttributeError:
            return None

    def _fset(self, obj, new_value):
        if obj.__getattribute__(self.public_name) is not None:
            obj.__getattribute__(self.public_name).__setattr__(self.container_name, None)

        try:
            if new_value.__getattribute__(self.container_name) is not None:
                warnings.warn('Controller already used by another `simulation`. Copied. Access from `simulation` object only',
                              Warning)
                new_value = copy.deepcopy(new_value)
        except AttributeError:
            pass

        try:
            new_value.__setattr__(self.container_name, obj)
        except AttributeError:
            pass

        obj.__setattr__(self.private_name, new_value)

    @property
    def private_name(self):
        return "_{}".format(self.f.__name__)

    @property
    def public_name(self):
        return self.f.__name__

    @property
    def container_name(self):
        return 'simulation'
