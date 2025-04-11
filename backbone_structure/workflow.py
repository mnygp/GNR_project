import taskblaster as tb


@tb.workflow
class Workflow:
    max_width = tb.var()

    @tb.dynamical_workflow_generator({'results': '*/*'})
    def generate_wfs(self):
        return tb.node('generate_wfs_from_list', input=self.max_width)

    @tb.task
    def collect_results(self):
        # Collect all the return_dict tasks from the results of the generated sub-workflows
        results = {}
        for (sub_workflow_name, sub_workflow_result in
             self.generate_wfs.results.items()):
                results[sub_workflow_name] = sub_workflow_result
        return results

    @tb.task
    def write_csv_task(self):
        # Collect all results from the workflow
        results = self.collect_results()

        # Process the results and write them to CSV
        return tb.node('write_results_to_csv', results_dict=results)


def workflow(runner):
    runner.run_workflow(Workflow(max_width=20))
