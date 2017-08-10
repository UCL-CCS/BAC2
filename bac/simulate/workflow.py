import subprocess
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
            sim = next(sim for sim in self._simulations if sim.is_ready)

            subprocess.run(sim.executable)

        print('Executing on {}'.format(self.resource))

    def preprocess_simulations(self):
        for *ensembles, simulation in product(*self.ensembles, self.simulations):

            sim = copy.deepcopy(simulation)
            self._simulations.append(sim)

            for ensemble in ensembles:
                ensemble.modifier(sim)

            prefix = Path(*(ens.name for ens in ensembles))
            self.path.joinpath(prefix).mkdir(parents=True, exist_ok=True)

            sim.restructure_paths_with_prefix(prefix=prefix)

            Encoder.encode(sim, self.path)

            p = subprocess.run(sim.preprocess_executable, shell=True, stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE, cwd=self.path)
            print(sim.preprocess_executable)
            print(p.stdout)
            if p.stderr: print(p.stderr)
            print()

    def __len__(self):
        return sum(1 if not x.is_finished else 0 for x in self._simulations)





