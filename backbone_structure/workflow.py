import taskblaster as tb


@tb.workflow
class Workflow:
    max_width = tb.var()

    @tb.dynamical_workflow_generator({'results': '*/*',
                                      'results_gap': '*/calculate_gap_task'})
    def generate_wfs(self):
        return tb.node('generate_wfs_from_list', input=self.max_width)

    @tb.task
    def write_csv_task(self):
        return tb.node('write_results_to_csv',
                       results_gap=self.generate_wfs.results_gap)


def workflow(runner):
    runner.run_workflow(Workflow(max_width=20))
