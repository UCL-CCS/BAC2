### A complex simulation example

First import the required modules
```python
import bac.simulate.gromacs as gmx
```

Instantiate the main `Run` class
```python
md_run = gmx.Run(integrator='md')
md_run.delta_time = 0.002
md_run.number_of_steps = 500000
```

Add constraints
---

```python
cc = gmx.ConstraintController()
cc.continuation = True
cc.type = gmx.ConstraintType.all_bonds
cc.algorithm = gmx.ConstraintAlgorithmType.lincs
md_run.constraint_controller = cc
```

Non bonded interactions
---

```python
nbc = gmx.NonBondedController()
nbc.coulomb_type = 'PME'
nbc.coulomb_cutoff = 1.0
md_run.non_bonded_controller = nbc
```
Temperature coupling
---

```python
tc = gmx.TemperatureController()
tc.type = 'V-rescale'
tc.time = 0.1
tc.temperature = 300
md_run.temperature_controller = tc
```

Pressure coupling
---

```python
pc = gmx.PressureController()
pc.type = 'Parrinello-Rahman'
pc.compressibility = 4.5e-5
pc.pressure = 1.0
pc.time = 2.0
md_run.pressure_controller = pc
```