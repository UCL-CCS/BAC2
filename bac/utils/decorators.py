from pathlib import Path
from enum import Enum
import copy

from bac.utils.pdb import PDBColumn
import warnings


class advanced_property(property):
    def __init__(self, *args, **kwargs):

        self._default = kwargs.get('default')
        self.type = kwargs.get('type')
        self.validator = kwargs.get('validator')
        self.warning_only = kwargs.get('warn', False)

        self.f = args[0] if args else None

        super(advanced_property, self).__init__(fget=self._fget, fset=self._fset)

    def __call__(self, f):
        self.f = f
        return self

    def _fget(self, obj):
        try:
            v = obj.__getattribute__(self.private_name)
        except AttributeError:
            v = None
        try:
            return v if v is not None else self.default(obj)
        except (AttributeError, TypeError):
            return None

    def _fset(self, obj, value):
        if value is None:
            obj.__setattr__(self.private_name, None)
        elif isinstance(value, self.type):
            value = self.type[0](value)
            self.validate(obj, value)
            obj.__setattr__(self.private_name, value)
        else:
            raise TypeError("{} must be of type {} instead of {}".format(self.public_name,
                                                                         ' or '.join(t.__name__ for t in self.type),
                                                                         type(value).__name__))

    @property
    def private_name(self):
        return "_{}".format(self.f.__name__)

    @property
    def public_name(self):
        return self.f.__name__

    def default(self, obj):
        if self._default is None: return None

        default_to_return = self._default(obj) if callable(self._default) else self._default

        return self.type[0](default_to_return)

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
            else:
                self._type = (value, )

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
