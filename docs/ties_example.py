import bac.simulate.namd as nd

from bac.simulate.coding import Encoder, Engine

tr = nd.Run()

tr.parameter_type_CHARMM = False
tr.amber = True
tr.parameters = '../build/complex.top'
tr.read_exclusions = False

tr.non_bonded_controller.excluded_interactions = 'scaled1-4'
tr.non_bonded_controller.one_four_scaling = 0.833333
tr.non_bonded_controller.cutoff = 12
tr.non_bonded_controller.switching = True
tr.non_bonded_controller.switching_distance = 10
tr.non_bonded_controller.pairlist_distance = 13.5
tr.non_bonded_controller.non_bonded_frequency = 1
tr.non_bonded_controller.full_elect_frequency = 2
tr.non_bonded_controller.steps_per_cycle = 10

tr.boundary_condition.wrap_water = True
tr.boundary_condition.wrap_all = True

tr.timestep = 2.0
tr.coordinates = '../build/complex.pdb'
tr.constraints = True
tr.temperature = 300

tr.temperature_controller.langevin = nd.LangevinDynamics(damping=5, temperature=200, hydrogen=False)

tr.binary_output = False
tr.output_name = '../equilibration/eq0'

print(Encoder.encode(tr, engine=Engine.namd))


if __name__ == '__main__':
    pass

