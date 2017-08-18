import bac.simulate.gromacs as gmx
from bac.simulate import Workflow, lambda_window, replica


class TiesGromacsEq(gmx.Run):
    def __init__(self):
        super(TiesGromacsEq, self).__init__()
        # Add further customizations...


class TiesGromacsMD(gmx.Run):
    def __init__(self):
        super(TiesGromacsMD, self).__init__()
        # Add further customizations...


class TiesWorkflow(Workflow):
    def __init__(self, resource, name):
        super(TiesWorkflow, self).__init__(resource, name)

        self.ensembles = [replica(5), lambda_window(13)]

        eq1 = TiesGromacsEq()
        eq2 = TiesGromacsEq()
        md = TiesGromacsMD()

        eq2.add_input_dependency(eq1)
        md.add_input_dependency(eq2)

        self.simulations = [eq1, eq2]

if __name__ == '__main__':
    # One can then do this:
    # Files to top, coord etc. needs to be set, not shown here...
    import bac.protocols.ties as ti
    wf = ti.TiesWorkflow(resource='lrz', name='test_run')
    wf.execute()
