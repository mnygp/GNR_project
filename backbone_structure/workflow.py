import taskblaster as tb


@tb.workflow
class Workflow:
    max_width = tb.var()

    @tb.dynamical_workflow_generator({'results': '*/*'})
    def generate_wfs(self):
        return tb.node('generate_wfs_from_list', input=self.max_width)


def workflow(runner):
    runner.run_workflow(Workflow(max_width=20))
