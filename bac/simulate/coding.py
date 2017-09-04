from enum import Enum
import yaml

from .basesimulation import BaseSimulation, Engine



class Encoder:
    """Encoding functionality. Singleton-like.

    """

    # These are simple lambda functions that turn
    # python types into specification correct
    # format. A notable example is the yes/no and
    # on/off bool problem.
    _encoders = {'yes_no': lambda x: ('no', 'yes')[x],
                 'on_off': lambda x: ('off', 'on')[x],
                 'enum': lambda x: x.value,
                 'iter': lambda x: ' '.join(str(xs) for xs in x),
                 'two_digit_iter': lambda xs: ' '.join(f'{x:.2f}'for x in xs)}

    _column_size = 25

    @classmethod
    def encode(cls, run: BaseSimulation, path):
        """

        Parameters
        ----------
        run
            Simulation class
        path: Path
            Location of configuration file

        """
        path = path if path.suffix == run.configuration_file_suffix else path / run.name.with_suffix(run.configuration_file_suffix)

        def serialize(obj, attributes):
            serial = ''
            for external_name, internals in attributes.items():
                # As a simplification, if the internal and external name are the same
                # one can omit the internal name in the YAML file. This results in
                # `internals` to be `None`. We can set `internals` to be `external_name`.
                if internals is None:
                    internals = external_name

                # Simple key value pair
                if isinstance(internals, str):
                    value = obj.__getattribute__(internals)
                    if value is not None:
                        serial += "{}{} = {}\n".format(external_name,' '*(cls._column_size-len(external_name)), value)

                # The value is a list containing the name and encoding mechanism
                # The list has to be [name, encoding] format. Nothing else works.
                elif isinstance(internals, list):
                    value = obj.__getattribute__(internals[0])
                    if value is not None:
                        encoder = cls._encoders[internals[1]]
                        serial += "{}{} = {}\n".format(external_name, ' '*(cls._column_size-len(external_name)),
                                                       encoder(value))

                # The value is a dictionary. Meaning that we have to go one level
                # down in the hierarchy. This function is then called recursively on
                # the object.
                elif isinstance(internals, dict):
                    sub_obj = obj.__getattribute__(external_name)

                    name_token = internals.pop('name_token', None)
                    if name_token is not None:
                        encoder = cls._encoders[name_token[1]]
                        serial += "{}{} = {}\n".format(name_token[0], ' '*(cls._column_size-len(name_token[0])),
                                                       encoder(sub_obj is not None))

                    if sub_obj is not None:
                        serial += serialize(sub_obj, internals)
            return serial+'\n'

        with open(path, mode='w') as conf:
            conf.write(serialize(run, cls._schema(Engine(run.engine_type.value))))

    # FIXME: setuptools for resource (yaml) allocation.
    @staticmethod
    def _schema(engine: Engine):
        if engine is Engine.namd:
            return yaml.load(open('/Users/kristofarkas/Developer/BAC2/bac/simulate/namd/namd_schema.yaml'))
        elif engine is Engine.gromacs:
            return yaml.load(open('/Users/kristofarkas/Developer/BAC2/bac/simulate/gromacs/gromacs_schema.yaml'))
