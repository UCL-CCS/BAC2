from pathlib import Path
import copy
import warnings
from typing import Callable, Union
from functools import partial

import numpy as np

from bac.utils.pdb import PDBColumn
from .versioned import Versioned


class advanced_property(property, Versioned):
    """ A `property` subclass with bells and whistles.

    This is a rather specific (yet common) implementation of the property mechanism. The getter returns the
    _<property name> or the default if that is None. The setter sets _<property name> to the value after
    validating it. If the value is set to None, then the getter *will* return the default. **In summary** you
    do not have to do any of the implementation! It **is** enough to ```def property_name(self): pass```.


    Examples
    --------


    """
    def __init__(self, *args, **kwargs):
        """

        Parameters
        ----------
        default
            The default value for the property. Can be a lambda function too. Self is passed at the only
            parameter of the lambda function to evaluate more complex default values.
        type: Union[type, Tuple(type)]
            This is a type or any function that takes anything supported and converts in into the desired
            correct type, or fails if the type or the value is not supported or correct.
        validator: lambda container_object, value -> bool
            This is a function that return a Bool representing whether the value is valid. Note that the object
            is passed to the function too, so one can make more complicated validation based on the state of the
            object.
        warn : bool
            If the validation fails, should it cause a warning or raise an Exception.
            Default: False
        available: StrictVersion
            The minimum version number that this property is supported on.

        Examples
        --------

        @advanced_property(type=(float, int), default=0.01, validator=lambda x, _: x >= 0, warn=True, available('5.0')
        def velocity(self): pass

        The property called `velocity` can be set as an float or int. Because float is the first in the list of types
        even is `self.velocity=2` is run `self.velocity` will return 2.0. Validator just checks non-negativeness. The
        system will only warn the user about a negative values, but will still set the property.


        """

        self._default = kwargs.get('default')
        self.type: Union[Callable, type] = kwargs.get('type')
        self.validator: Callable = kwargs.get('validator')
        self.warning_only = kwargs.get('warn', False)
        self.warning_message = kwargs.get('warning_message')
        self.version = kwargs.get('available')
        self.post_set_processor = None

        self.f = args[0] if args else None

        super(advanced_property, self).__init__(fget=self._fget, fset=self._fset)

    def __call__(self, f):
        self.f = f
        return self

    def _fget(self, obj):

        self.version_check(obj)

        try:
            v = obj.__getattribute__(self.private_name)
        except AttributeError:
            v = None

        try:
            return v if v is not None else self.default(obj)
        except (AttributeError, TypeError) as e:
            warnings.warn(f'{e}; {self.f.__name__}', Warning)
            return None

    def _fset(self, obj, value):
        self.version_check(obj)

        if value is None:
            obj.__setattr__(self.private_name, None)
        else:
            value = value if isinstance(self.type, type) and isinstance(value, self.type) else self.type(value)
            try:
                self.validate(obj, value)
            except AttributeError as e:
                warnings.warn(f'{e}')
            obj.__setattr__(self.private_name, value)

        if self.post_set_processor:
            self.post_set_processor(obj, value)

    @property
    def private_name(self):
        return "_{}".format(self.f.__name__)

    @property
    def public_name(self):
        return self.f.__name__

    def default(self, obj):
        if self._default is None:
            return None

        default_to_return = self._default(obj) if callable(self._default) else self._default

        if default_to_return is None:
            return None
        elif isinstance(self.type, type) and isinstance(default_to_return, self.type):
            return default_to_return
        else:
            return self.type(default_to_return)

    def validate(self, obj, value):
        if self.validator is not None and self.validator(obj, value) is False:
            message = self.warning_message or "Setting {} to {} does not fulfill restrictions.".format(self.public_name, value)
            if self.warning_only:
                warnings.warn(message, Warning)
            else:
                raise ValueError(message)

    def version_check(self, obj):
        if isinstance(obj, Versioned) and self.version is not None and obj.version is not None:
            if self.version > obj.version:
                warnings.warn('This is not supported on current version!', Warning)

    @property
    def f(self):
        return self._f

    @f.setter
    def f(self, f):
        self.__doc__ = f.__doc__
        self._f = f

    def post_set_processing(self, f):
        self.post_set_processor = f
        return self


class pathlike(advanced_property):
    def __init__(self, *args, **kwargs):
        super(pathlike, self).__init__(type=Path, *args, **kwargs)


class pdbcolumn(advanced_property):
    def __init__(self, *args, **kwargs):
        kwargs['type'] = PDBColumn
        if kwargs.get('default') is None:
            kwargs['default'] = 'O'
        super(pdbcolumn, self).__init__(*args, **kwargs)


class decimal(advanced_property):
    def __init__(self, *args, **kwargs):
        super(decimal, self).__init__(type=np.float, *args, **kwargs)


class positive_decimal(decimal):
    def __init__(self, *args, **kwargs):
        old_validator = kwargs.get('validator', lambda o, v: True)
        kwargs['validator'] = lambda o, v: old_validator(o, v) and v >= 0
        super(positive_decimal, self).__init__(*args, **kwargs)


class integer(advanced_property):
    def __init__(self, *args, **kwargs):
        super(integer, self).__init__(type=np.int, *args, **kwargs)


class positive_integer(integer):
    def __init__(self, *args, **kwargs):
        old_validator = kwargs.get('validator', lambda o, v: True)
        kwargs['validator'] = lambda o, v: old_validator(o, v) and v > 0
        super(positive_integer, self).__init__(*args, **kwargs)


class non_negative_integer(integer):
    def __init__(self, *args, **kwargs):
        old_validator = kwargs.get('validator', lambda o, v: True)
        kwargs['validator'] = lambda o, v: old_validator(o, v) and v >= 0
        super(non_negative_integer, self).__init__(*args, **kwargs)


class boolean(advanced_property):
    def __init__(self, *args, **kwargs):
        super(boolean, self).__init__(type=bool, *args, **kwargs)


class float_vector(advanced_property):
    def __init__(self, *args, **kwargs):
        super(float_vector, self).__init__(type=partial(np.array, dtype=np.float, ndmin=1), *args, **kwargs)


class back_referenced(property):
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
