from enum import Enum
import yaml


class Engine(Enum):
    namd = 'namd'
    gromacs = 'gromacs'


class Encoder:

    _encoders = {'yes_no': lambda x: ('no', 'yes')[x],
                 'on_off': lambda x: ('off', 'on')[x],
                 'enum': lambda x: x.name}

    @classmethod
    def encode(cls, run, engine):

        def serialize(obj, attributes):
            s = ""
            for external_name, internals in attributes.items():
                if isinstance(internals, str):
                    value = obj.__getattribute__(internals)
                    if value is not None:
                        s += "{} {}\n".format(external_name, value)
                elif isinstance(internals, list):
                    value = obj.__getattribute__(internals[0])
                    if value is not None:
                        encoder = cls._encoders[internals[1]]
                        s += "{} {}\n".format(external_name, encoder(value))
                elif isinstance(internals, dict):
                    sub_obj = obj.__getattribute__(external_name)

                    name_token = internals.pop('name_token', None)
                    if name_token is not None:
                        encoder = cls._encoders[name_token[1]]
                        s += "{} {}\n".format(name_token[0], encoder(sub_obj is not None))

                    if sub_obj is not None:
                        s += serialize(sub_obj, internals)
            return s

        return serialize(run, cls._schema(engine))

    @staticmethod
    def _schema(engine):
        if engine is Engine.namd:
            return yaml.load(open('/Users/kristofarkas/Developer/BAC2/bac/simulate/namd/namd_schema.yaml'))
        elif engine is Engine.gromacs:
            return yaml.load(open('/Users/kristofarkas/Developer/BAC2/bac/simulate/gromacs/gromacs_schema.yaml'))
