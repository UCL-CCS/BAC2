# The *simulate* modules

BAC2.0 support Gromacs and NAMD simulation engines.
While they have their separate class structure, you can use them 
interchangeably in most cases. You can also convert between the two.

### Setup

All md abstractions have as their main class a `run_controller.py` module.
Whenever you want to run a new simulation you have to instantiate a new `Run()`
class.  
If you are creating a workflow, you most likely will *__only__* need to `deepcopy` 
the first instance, and modify some parameter and the input/output.


**Run** _main class_   
- TemperatureController  
- PressureController  
- NonBondedController  
- ConstraintController

```python
import copy
import bac.simulate.namd as namd

annealing = namd.Run(temperature=300.0)
annealing.temperature_controller.langevin = namd.LangevinDynamics(damping=0.4)

annealing.coordinates = "../bace1.pdb"
# OR
annealing.input = bac.build.SOMECLASS()

equilibrate = copy.deepcopy(annealing)
equilibrate.temperature = 310.0

# This sets all the input file paths to the output paths of the
# pervious md simulation.
equilibrate.input = annealing


operation_queue = OperationQueue(resource='archer')

equilibrate.add_dependency(annealing)

operation_queue.add_operation(equilibrate) 
operation_queue.add_operation(annealing)

# Asynchronous execution 
operation_queue.execute()
```

