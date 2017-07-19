from enum import Enum

from bac.utils.decorators import advanced_property, integer, decimal

class PressureCouplingType(Enum):
    no = 'no'
    berendsen = 'berendsen'
    parrinello_rahman = 'Parrinello-Rahman'
    mttk = 'MTTK'

class IsotropyType(Enum):
    isotropic = 'isotropic'
    semi_isotropic = 'semiisotropic'
    anisotropic = 'anisotropic'
    surface_tension = 'surface_tension'


class PressureController:

    def __init__(self, **kwargs):
        self.type = kwargs.get('type')
        self.isotropy = kwargs.get('isotropy')
        self.frequency = kwargs.get('frequency')
        self.time = kwargs.get('time')
        self.compressibility = kwargs.get('compressibility')
        self.pressure = kwargs.get('pressure')
        self.reference_coordinate_scaling = kwargs.get('reference_coordinate_scaling')

    @advanced_property(type=PressureCouplingType, default=PressureCouplingType.no)
    def type(self): pass

    @advanced_property(type=IsotropyType, default=IsotropyType.isotropic)
    def isotropy(self): pass


