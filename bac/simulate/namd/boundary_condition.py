from collections import namedtuple

from bac.utils.decorators import *


Point3D = namedtuple('Point3D', 'x y z')


class BoundaryCondition:

    @property
    def is_periodic(self):
        return self is PeriodicBoundaryCondition


class PeriodicBoundaryCondition(BoundaryCondition):

    def __init__(self, **kwargs):
        self.cell_basis_vector_1 = kwargs.get('cell_basis_vector_1')
        self.cell_basis_vector_2 = kwargs.get('cell_basis_vector_2')
        self.cell_basis_vector_3 = kwargs.get('cell_basis_vector_3')
        self.cell_origin = kwargs.get('cell_origin')

        self.extended_system = kwargs.get('extended_system')

        self.XST_file = kwargs.get('XST_file')
        self.XST_frequency = kwargs.get('XST_frequency')

        self.wrap_water = kwargs.get('wrap_water')
        self.wrap_all = kwargs.get('wrap_all')
        self.wrap_nearest = kwargs.get('wrap_nearest')

    @advanced_property(type=Point3D, default=(0, 0, 0))
    def cell_basis_vector_1(self): pass

    @advanced_property(type=Point3D, default=(0, 0, 0))
    def cell_basis_vector_2(self): pass

    @advanced_property(type=Point3D, default=(0, 0, 0))
    def cell_basis_vector_3(self): pass

    @advanced_property(type=Point3D, default=(0, 0, 0))
    def cell_origin(self): pass

    @file
    def extended_system(self): pass

    @file(default=lambda s: s.run.output_name.with_suffix('.xst'))
    def XST_file(self): pass

    @positive_integer
    def XST_frequency(self): pass

    @boolean(default=False)
    def wrap_water(self):pass

    @boolean(default=False)
    def wrap_all(self): pass

    @boolean(default=False)
    def wrap_nearest(self): pass

# TODO There are two other types of boundary conditions. Please implement.
