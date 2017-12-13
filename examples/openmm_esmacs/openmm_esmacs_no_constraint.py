import simtk.openmm.app as app
import simtk.openmm as mm
import simtk.unit as u

import numpy as np

import time

# System

prmtop = app.AmberPrmtopFile('complex.top')
pdb = app.PDBFile('complex.pdb')
system = prmtop.createSystem(nonbondedMethod=app.PME,
                             constraints=app.HBonds,
                             nonbondedCutoff=12*u.angstroms,
                             switchDistance=10*u.angstroms)

total_steps = 100000 # Reducing for testing purposes from 3M

# Integrator, forces

integrator = mm.LangevinIntegrator(50*u.kelvin, 5/u.picosecond, 0.002*u.picoseconds)

barostat = mm.MonteCarloBarostat(1.0*u.bar, 300*u.kelvin)

# Initialise simulation and positions

simulation = app.Simulation(prmtop.topology, system, integrator)

simulation.context.setPositions(pdb.positions)
simulation.context.setVelocitiesToTemperature(50*u.kelvin)

# Minimise

print('Minimising...')

initial_time = time.time()

simulation.minimizeEnergy(maxIterations=1000)

# Heating

print('Heating...')

for temperature in np.linspace(50, 300, 251)*u.kelvin:
    integrator.setTemperature(temperature)
    simulation.step(100)
    print(temperature, end=' ')

simulation.step(5000)

final_time = time.time()
elapsed_time = (final_time - initial_time) * u.seconds
print('Completed %8d steps in %8.3f s' % (total_steps, elapsed_time / u.seconds))
print('Potential energy is %.3f kcal/mol' % (simulation.context.getState(getEnergy=True).getPotentialEnergy() / u.kilocalories_per_mole))

# Equilibrate

print('Benchmarking...')
print('Equilibrate...')
initial_time = time.time()

simulation.system.addForce(barostat)
simulation.step(int(total_steps*1/3))
simulation.saveState('equilibrated.xml')

# Production
print('Production...')
simulation.reporters.append(app.CheckpointReporter('checkpnt.chk', int(total_steps/10)))
simulation.reporters.append(app.DCDReporter('simulation.dcd', int(total_steps/400)))
simulation.reporters.append(app.StateDataReporter('simulation.log', int(total_steps/200), step=True, time=True, totalEnergy=True))
simulation.step(int(total_steps*2/3))

final_time = time.time()
elapsed_time = (final_time - initial_time) * u.seconds
simulated_time = total_steps * integrator.getStepSize()
performance = (simulated_time / elapsed_time)
print('Completed %8d steps in %8.3f s : performance is %8.3f ns/day' % (total_steps, elapsed_time / u.seconds, performance / (u.nanoseconds/u.day)))
print('Final potential energy is %.3f kcal/mol' % (simulation.context.getState(getEnergy=True).getPotentialEnergy() / u.kilocalories_per_mole))
