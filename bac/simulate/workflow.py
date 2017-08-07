import subprocess
import uuid
import copy
from pathlib import Path
from itertools import product
import random
from bac.simulate.coding import Encoder


class Workflow:

    def __init__(self, resource, name):
        self.resource = resource
        self.path = Path(name + '_' + str(random.randint(1000, 9999)))

        self.simulations = []
        self._simulations = []
        self.ensembles = []

    def add_simulation(self, simulation):
        self.simulations.append(simulation)

    def execute(self):

        self.preprocess_simulations()

        while len(self):
            op = next(op for op in self._simulations if op.is_ready)

            subprocess.run(op.execute_path)

        print('Executing on {}'.format(self.resource))

    def preprocess_simulations(self):
        for *ensembles, simulation in product(*self.ensembles, self.simulations):

            sim = copy.deepcopy(simulation)
            self._simulations.append(sim)

            for ensemble in ensembles:
                ensemble.modifier(sim)

            full_path = self.path.joinpath(*(ens.name for ens in ensembles))
            full_path.mkdir(parents=True, exist_ok=True)

            sim.restructure_paths_with_prefix(prefix=full_path)

            Encoder.encode(sim, full_path)

            p = subprocess.run(sim.preprocess_executable, shell=True, stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)
            # print(sim.preprocess_executable)
            # print(p.stdout)
            # if p.stderr: print(p.stderr)
            # print()

    def __len__(self):
        return sum(1 if not x.is_finished else 0 for x in self._simulations)





