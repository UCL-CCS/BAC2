import subprocess
import uuid
from pathlib import Path
from itertools import product
import copy
from bac.simulate.coding import Encoder


class Workflow:

    def __init__(self, resource, name):
        self.resource = resource
        self.name = name or str(uuid.uuid4())
        self.path = Path(self.name)
        self.simulations = []
        self._simulations = []
        self.modifiers = []

    def add_simulation(self, simulation):
        self.simulations.append(simulation)

    def execute(self):

        self.preprocess_simulations()

        while len(self):
            op = next(op for op in self._simulations if op.is_ready)

            subprocess.run(op.execute_path)

        print('Executing on {}'.format(self.resource))

    def preprocess_simulations(self):
        for mods in product(*self.modifiers):

            full_path = self.path
            sims = copy.deepcopy(self.simulations)

            for mod in mods:
                full_path = full_path / mod[1]
                for sim in sims:
                    mod[0](sim)

            full_path.mkdir(parents=True)
            for sim in sims:
                Encoder.encode(sim, full_path)

                for p in sim.executable_path:
                    print(' '.join(str(a) for a in p))

            self._simulations += sims

    def __len__(self):
        return sum(1 if not x.is_finished else 0 for x in self._simulations)





