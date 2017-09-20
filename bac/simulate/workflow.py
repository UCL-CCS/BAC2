import subprocess
import copy
from pathlib import Path
from itertools import product
import random
from bac.simulate.coding import Encoder
from bac.simulate.basesimulation import BaseSimulation
from bac.simulate.ensemble import BaseEnsembleIterator

from typing import List


class Workflow:
    """

    Attributes
    ----------
    simulations: list
        The list of simulations in the workflow. This is generic and not dependent
        on the ensembles. The actual simulation list, which is for all ensembles
        is generated at `execute` time.
    ensembles: list
        A list of ensemble mechanism used to duplicated the simulations.
    """

    def __init__(self, resource, r_dir: Path):
        """

        Parameters
        ----------
        resource: str
            The supercomputer the workflow will run on.
        dir: path like
            The name of the simulation. This will be the master folder name too.

        Methods
        -------


        """

        r_dir = Path(r_dir)

        if r_dir.exists() and r_dir.is_dir():
            import shutil
            shutil.rmtree(r_dir)

        self.resource: str = resource
        self.path: Path = r_dir
        self.simulations: List[BaseSimulation] = []
        self._simulations: List[BaseSimulation] = []
        self.ensembles: List[BaseEnsembleIterator] = []

    def add_simulation(self, simulation: BaseSimulation):

        """

        Parameters
        ----------
        simulation: BaseSimulation

        Returns
        -------

        """
        self.simulations.append(simulation)

    def execute(self):

        self.preprocess_simulations()

        while len(self):
            sim = next(sim for sim in self._simulations if sim.is_ready)

            subprocess.run(sim.executable)

        print('Executing on {}'.format(self.resource))

    def preprocess_simulations(self):
        """Run pre-processing tasks for the simulations.

        Parameters
        ----------
        execute: bool
            Execute the pre processing step on the shell. If `False` then the
            executable is printed to stdout.

        """
        for *ensembles, simulation in product(*self.ensembles, self.simulations):

            sim = copy.deepcopy(simulation)
            self._simulations.append(sim)

            for ensemble in ensembles:
                ensemble.fn(sim)

            prefix = Path(*(ens.path_name for ens in ensembles))
            self.path.joinpath(prefix).mkdir(parents=True, exist_ok=True)

            sim.restructure_paths_with_prefix(prefix=prefix)

            Encoder.encode(sim, self.path)

        self.write_generic_bash_executable(path=self.path)

    def write_generic_bash_executable(self, path: Path):

        if not path.is_dir():
            raise NotADirectoryError

        with open(path / 'workflow.sh', mode='w') as wf, open(path / 'run.sh', mode='w') as rn:

            wf.write('#!/usr/bin/env bash\n\n')
            rn.write('#!/usr/bin/env bash\n\n'
                     '#PBS -l nodes=%%%:ppn=32:xe\n'
                     '#PBS -l walltime=%%:%%:00\n'
                     'module swap PrgEnv-cray PrgEnv-gnu\n'
                     'export OMP_NUM_THREADS=1\n'
                     'cd $PBS_O_WORKDIR\n')

            for index, ens in enumerate(self.ensembles, start=1):
                wf.write(f"{ens.name.upper()}=${index}\n")

            wf.write('\n')

            for simulation in self.simulations:
                sim = copy.deepcopy(simulation)

                prefix = Path(*(ens.generic_path_name for ens in self.ensembles))

                sim.restructure_paths_with_prefix(prefix=prefix)

                wf.write(sim.preprocess_executable+'\n')
                wf.write(sim.executable+'\n\n')

            for ens in product(*self.ensembles):
                rn.write(f"bash workflow.sh {' '.join(str(e.iterator_state) for e in ens)} &\n")

            rn.write('wait\nexit 0\n')

    def __len__(self):
        return sum(1 if not x.is_finished else 0 for x in self._simulations)



