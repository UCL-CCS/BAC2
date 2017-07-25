# The *simulate* modules

BAC2.0 support Gromacs and NAMD simulation engines.
While they have their separate class structure, you can use them 
interchangeably in most cases. You can also convert between the two.

### Setup

All md abstractions have as their main class a `Run` object inside the `run_controller` module.
Whenever you want to run a new simulation you have to instantiate a new `Run`
class.  
If you are creating a workflow, you most likely will *__only__* need to `deepcopy` 
the first instance, and modify some parameter and the input/output. *__Note__: I think an even 
 better way to handle reference semantics problem is by copying the passed Run class when adding 
 it to the* `operation_queue`.


**Run** – the _**main** class_   
- TemperatureController  
- PressureController  
- NonBondedController  
- ConstraintController

```python
import copy
import bac.simulate.namd as namd
from bac.simulate.operation_queue import OperationQueue
from bac.simulate.operation import Operation

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


ties_bace1 = OperationQueue(resource='archer')

for _ in range(5):
    op_annealing = Operation(run=annealing)
    
    op_equilibrate = Operation(run=equilibrate)
    op_equilibrate.add_dependency(op_annealing)
    
    ties_bace1.add_operation(op_annealing, op_equilibrate)

# Asynchronous execution
ties_bace1.execute()
```
end
