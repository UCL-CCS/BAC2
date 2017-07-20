#The *simulate* modules

BAC2.0 support Gromacs and NAMD simulation engines.
While they have their separate class structure, you can use them 
interchangeably in most cases. You can also convert between the two.

###Setup

All md abstractions have as their main class a `run_controller.py` module.
Whenever you want to run a new simulation you have to instantiate a new `Run()`
class.  
If you are creating a workflow, you most likely will *only* need to `deepcopy` 
the first instance, and modify some parameter and the input/output.


**Run** _main class_   
- TemperatureController  
- PressureController  
- NonBondedController  
- ConstraintController

```python
import copy
import bac.simulate.namd as nd

annealing = nd.Run(temperature=300.0)
annealing.coordinates = "../bace1.pdb"

equilibrate = copy.deepcopy(annealing)
equilibrate.temperature = 310.0

operation_queue = OperationQueue(resource='archer')

equilibrate.add_dependency(annealing)

# Order doesn't matter. Dependencies will be fullfilled
operation_queue.add_operation(equilibrate, annealing) 

# Asynchronous execution 
operation_queue.execute()
```

