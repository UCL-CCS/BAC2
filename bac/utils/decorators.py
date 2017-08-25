from pathlib import Path
from enum import Enum
import copy

from bac.utils.pdb import PDBColumn
import warnings

from .versioned import Versioned


class advanced_property(property, Versioned):
    """ A `property` subclass with bells and whistles.

    This is a rather specific (yet common) implementation of the property mechanism. The getter return the
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
            Possible types that property can be set to. **Important**: if more than one (1) type is set,
            then the first one will be the major, and this first one *must* be instantiatable from the
            other types.
        validator: lambda value, object -> bool
            This is a function that return a Bool representing wether the value is valid. Note that the object
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
        self.type = kwargs.get('type')
        self.validator = kwargs.get('validator')
        self.warning_only = kwargs.get('warn', False)
        self.version = kwargs.get('available')

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
        except (AttributeError, TypeError):
            return None

    def _fset(self, obj, value):
        self.version_check(obj)

        if value is None:
            obj.__setattr__(self.private_name, None)
        elif isinstance(value, self.type):
            value = self.convert_to_default_type(value)
            self.validate(obj, value)
            obj.__setattr__(self.private_name, value)
        else:
            raise TypeError("{} must be of type {} NOT {}".format(self.public_name, ' or '.join(t.__name__ for t in self.type), type(value).__name__))

    @property
    def private_name(self):
        return "_{}".format(self.f.__name__)

    @property
    def public_name(self):
        return self.f.__name__

    def default(self, obj):
        if self._default is None: return None

        default_to_return = self._default(obj) if callable(self._default) else self._default

        return self.convert_to_default_type(default_to_return)

    @property
    def type(self):
        return self._type

    @type.setter
    def type(self, value):
        if isinstance(value, tuple):
            self._type = value
        else:
            if issubclass(value, Enum):
                self._type = (value, str)
            elif value is list:
                self._type = (value, str, int, float)
            else:
                self._type = (value, )

    def convert_to_default_type(self, value):
        default_type = self.type[0]
        if default_type is list and not isinstance(value, list):
            value = [value]

        return default_type(value)

    def validate(self, obj, value):
        if self.validator is not None and self.validator(value, obj) is False:
            message = "Setting {} to {} does not fulfill restrictions.".format(self.public_name, value)
            if self.warning_only:
                warnings.warn(message, Warning)
            else:
                raise ValueError(message)

    @property
    def validator(self):
        return self._validator

    @validator.setter
    def validator(self, value):
        if value is None:
            self._validator = None
            return

        new_validator = value
        argument_count = value.__code__.co_argcount if value.__code__.co_argcount > 0 else None
        self._validator = lambda *a: new_validator(*a[:argument_count])

    def version_check(self, obj):
        if isinstance(obj, Versioned) and self.version is not None and obj.version is not None:
            if self.version > obj.version:
                warnings.warn('This is not supported on current version!', Warning)

class file(advanced_property):
    def __init__(self, *args, **kwargs):
        super(file, self).__init__(type=(Path, str), *args, **kwargs)


class column(advanced_property):
    def __init__(self, *args, **kwargs):
        super(column, self).__init__(type=PDBColumn, default=PDBColumn.O, *args, **kwargs)


class positive_decimal(advanced_property):
    def __init__(self, *args, **kwargs):
        old_validator = kwargs.get('validator', lambda *y: True)
        kwargs['validator'] = lambda *x: old_validator(*x) and x[0] >= 0
        super(positive_decimal, self).__init__(type=(float, int, str), *args, **kwargs)


class decimal(advanced_property):
    def __init__(self, *args, **kwargs):
        super(decimal, self).__init__(type=(float, int, str), *args, **kwargs)

class integer(advanced_property):
    def __init__(self, *args, **kwargs):
        super(integer, self).__init__(type=(int, str), *args, **kwargs)

class positive_integer(advanced_property):
    def __init__(self, *args, **kwargs):
        old_validator = kwargs.get('validator', lambda *y: True)
        kwargs['validator'] = lambda *x: old_validator(*x) and x[0] >= 0
        super(positive_integer, self).__init__(type=(int, str), *args, **kwargs)


class boolean(advanced_property):
    def __init__(self, *args, **kwargs):
        super(boolean, self).__init__(type=bool, *args, **kwargs)


class back_referenced(property):
    """Decorator for attributes that reference their owner.

    The decorated object has a property called `run` that references back to the
    main `Run` object. This is important for example when validators have to access
    properties from other parts of the object tree.

    Examples
    --------

    @back_referenced
    def pressure_controller(self): pass

    # Then in the class definition access self.run:

    class PressureController:

        @decimal(validator:lambda x, s: x < s.run.pressure)
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
        except:
            return None

    def _fset(self, obj, new_value):
        if obj.__getattribute__(self.public_name) is not None:
            obj.__getattribute__(self.public_name).__setattr__(self.container_name, None)

        try:
            if new_value.__getattribute__(self.container_name) is not None:
                warnings.warn('Controller already used by another `run`. Copied. Access from `run` object only',
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
        return 'run'
