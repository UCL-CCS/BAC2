from enum import Enum

from bac.utils.decorators import advanced_property, integer, decimal
from bac.simulate.gromacs.integrator import Integrator
from bac.simulate.coding import Encodable


class PressureCouplingType(Enum):
    no = 'no'
    berendsen = 'berendsen'
    parrinello_rahman = 'Parrinello-Rahman'
    mttk = 'MTTK'

    @classmethod
    def _missing_(cls, value):
        if value is False:
            return cls.no
        else:
            super()._missing_(value)


class IsotropyType(Enum):
    isotropic = 'isotropic'
    semi_isotropic = 'semiisotropic'
    anisotropic = 'anisotropic'
    surface_tension = 'surface_tension'


class ReferenceCoordinateScalingType(Enum):
    no = 'no'
    all = 'all'
    com = 'com'

    @classmethod
    def _missing_(cls, value):
        if value is False:
            return cls.no
        else:
            super()._missing_(value)


class PressureController(Encodable):

    def __init__(self, **kwargs):
        self.type = kwargs.get('type')
        self.isotropy = kwargs.get('isotropy')
        self.frequency = kwargs.get('frequency')
        self.time = kwargs.get('time')
        self.compressibility = kwargs.get('compressibility')
        self.pressure = kwargs.get('pressure')
        self.reference_coordinate_scaling = kwargs.get('reference_coordinate_scaling')

    def __eq__(self, other):
        if isinstance(other, PressureController):
            return (self.type == other.type and
                    self.isotropy == other.isotropy and
                    self.frequency == other.frequency and
                    self.time == other.time and
                    self.compressibility == other.compressibility and
                    self.pressure == other.pressure and
                    self.reference_coordinate_scaling == other.reference_coordinate_scaling)

        return False

    def __str__(self):
        return "{} {} pressure coupling\n" \
               "frequency {}, every {} ps\n" \
               "pressure: {} bar\n" \
               "compressibility: {} bar^-1".format(self.isotropy.name.capitalize(),
                                                   self.type.name.capitalize(),
                                                   self.frequency,
                                                   self.time,
                                                   self.pressure,
                                                   self.compressibility)

    def _repr_html_(self):
        return "<h1>{} {} pressure coupling</h1>\n" \
               "<p>Frequency {}, every {} ps</p>" \
               "<p>Pressure: {} bar</p>" \
               "<p>Compressibility: {} bar^-1</p>".format(self.isotropy.name.capitalize(),
                                                          self.type.name.capitalize(),
                                                          self.frequency,
                                                          self.time,
                                                          self.pressure,
                                                          self.compressibility)

    @advanced_property(type=PressureCouplingType, default=PressureCouplingType.no)
    def type(self): pass

    @advanced_property(type=IsotropyType, default=IsotropyType.isotropic)
    def isotropy(self): pass

    def _frequency_default(self):
        if self.simulation.integrator in (Integrator.md_vv, Integrator.md_vv_avek):
            return 1
        else:
            return self.simulation.non_bonded_controller.neighbor_list_update_frequency if self.simulation.non_bonded_controller.neighbor_list_update_frequency > 0 else 10

    @integer(default=_frequency_default)
    def frequency(self): pass

    @decimal(default=1)
    def time(self): pass

    def _group_count_validator(self, x):
        if self.isotropy is IsotropyType.isotropic:
            return len(x) == 1
        elif self.isotropy is IsotropyType.semi_isotropic:
            return len(x) == 2
        elif self.isotropy is IsotropyType.anisotropic:
            return len(x) == 6
        elif self.isotropy is IsotropyType.surface_tension:
            return len(x) == 2

    @advanced_property(type=list, default=[], validator=_group_count_validator)
    def compressibility(self): pass

    @advanced_property(type=list, default=[], validator=_group_count_validator)
    def pressure(self): pass

    @advanced_property(type=ReferenceCoordinateScalingType, default=ReferenceCoordinateScalingType.no)
    def reference_coordinate_scaling(self): pass
