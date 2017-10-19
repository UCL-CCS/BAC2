import yaml
from functools import reduce
from pathlib import Path

_encoders = {'yes_no': lambda x: ('no', 'yes')[x],
             'on_off': lambda x: ('off', 'on')[x],
             'enum': lambda x: x.value,
             'iter': lambda x: ' '.join(str(xs) for xs in x),
             'two_digit_iter': lambda xs: ' '.join(f'{x:.2f}'for x in xs)}


_column_size = 25


class Encodable:

    @staticmethod
    def _get_schema(class_name):
        return Path(__file__).parent.joinpath('schemas',f'{class_name[21:-2]}.yaml')

    def encode(self, path=None, suffix=None):

        # This is an ugly hack. So the yaml files are named as if you were to write out the
        # class' name, with the module hierarchy. But printing that out puts a < class ' in front
        # so we have to take that out, plus we also take out bac.simulate.
        # Also, this has to be set as a resource in the setup.py file, obv referencing it like this
        # is not cool.

        def update(d: dict, other: dict) -> dict: d.update(other); return d

        attributes = reduce(update, (yaml.load(open(self._get_schema(str(class_names))))
                                     for class_names in (self.__class__, *self.__class__.__bases__) if (issubclass(class_names, Encodable)) and class_names is not Encodable), {})

        s = self._serialize(attributes)

        if suffix:
            s += suffix

        if path is None:
            return s
        else:
            with open(path, mode='w') as f:
                f.write(s)

    def _serialize(self, attributes):
        serial = ''
        for external_name, internals in attributes.items():
            # As a simplification, if the internal and external name are the same
            # one can omit the internal name in the YAML file. This results in
            # `internals` to be `None`. We can set `internals` to be `external_name`.
            if internals is None:
                internals = [external_name]

            if isinstance(internals, str):
                internals = [internals]

            if not isinstance(internals, list):
                raise SyntaxError('YAML schema is invalid!')

            if internals[0] == 'active_indicator':
                value = True
            else:
                value = self.__getattribute__(internals[0])

            if isinstance(value, Encodable):
                serial += value.encode()
            elif value is not None:
                serial += self._format_line(external_name, value if len(internals) == 1 else _encoders[internals[1]](value))

        return serial + '\n'

    @staticmethod
    def _format_line(this: str, that: str, separator=' = ', align=True, terminator='\n'):

        this = str(this)

        if align:
            this += ' '*(_column_size - len(this))
        this += separator + str(that) + terminator

        return this




