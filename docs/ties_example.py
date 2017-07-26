import bac.simulate.namd as nd

md = nd.Run()

md.parameter_type_CHARMM = False
md.amber = True
md.parameters = '../build/complex.top'
md.read_exclusions = False

nbc = md.non_bonded_controller
nbc.excluded_interactions = 'scaled1-4'
nbc.one_four_scaling = 0.833333
nbc.cutoff = 12
nbc.switching = True
nbc.switching_distance = 10
nbc.pairlist_distance = 13.5
nbc.non_bonded_frequency = 1
nbc.full_elect_frequency = 2
nbc.steps_per_cycle = 10
nbc.pme = nd.PME(grid_spacing=1.0)

bc = md.boundary_condition
bc.cell_basis_vector_1 = (87.425, 0.000, 0.000)
bc.cell_basis_vector_2 = (0.000, 97.862, 0.000)
bc.cell_basis_vector_3 = (0.000, 0.000, 76.087)
bc.cell_origin = (-0.071, -0.225, -0.355)
bc.wrap_water = True
bc.wrap_all = True

cc = md.constraints
cc.bond_constraint = nd.BondConstraint(bonds='all', tolerance=0.00001, iterations=100)
cc.harmonic_constraint = nd.HarmonicConstraint(exponent=2, reference_position_file='../build/complex.pdb',
                                               force_constant_file='../constraint/f4.pdb', force_constant_column='B')

md.timestep = 2.0
md.coordinates = '../build/complex.pdb'
md.temperature = 300

md.temperature_controller.langevin = nd.LangevinDynamics(damping=5, temperature=300, hydrogen=False)

md.binary_output = False
md.output_name = '../equilibration/eq0'

fe = nd.FreeEnergyController(type='ti', file='../build/tags.pdb', column='B', decouple=True)

fe.output_name = '../equilibration/eq0.alch'
fe.output_frequency = 1000
fe.van_der_waals_shift_coefficient = 5
fe.electronic_interaction_start_lambda = 0.45
fe.van_der_waals_end_lambda = 1.0

md.free_energy_controller = fe

print(md)


if __name__ == '__main__':
    pass

