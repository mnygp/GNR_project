import taskblaster as tb


@tb.workflow
class Workflow:
    max_m = tb.var()

    @tb.dynamical_workflow_generator({'n_3_results':
                                      '*/*',
                                      'ldos_n_3_results':
                                      '*/ldos_edge_state_task'})
    def wfs_n_3(self):
        return tb.node('generate_wfs', m_max=self.max_m, n=3)

    @tb.dynamical_workflow_generator({'n_1_results':
                                      '*/*',
                                      'ldos_n_1_results':
                                      '*/ldos_edge_state_task'})
    def wfs_n_1(self):
        return tb.node('generate_wfs', m_max=self.max_m, n=1)


def workflow(runner):
    runner.run_workflow(Workflow(max_m=6))
