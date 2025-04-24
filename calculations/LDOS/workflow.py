import taskblaster as tb


@tb.workflow
class Workflow:
    max_repeat = tb.var()

    @tb.dynamical_workflow_generator({'results': '*/*',
                                      'pre_relax':
                                      '*/ldos_pre_relaxation_task',
                                      'post_relax':
                                      '*/ldos_post_relaxation_task'})
    def generate_wfs(self):
        return tb.node('generate_wfs_from_list', input=self.max_repeat)


def workflow(runner):
    runner.run_workflow(Workflow(max_repeat=5))
