import taskblaster as tb
from pathlib import Path
import csv


@tb.workflow
class Workflow:
    max_width = tb.var()

    @tb.dynamical_workflow_generator({'results': '*/*',
                                      'results_gap': '*/calculate_gap_task'})
    def generate_wfs(self):
        return tb.node('generate_wfs_from_list', input=self.max_width)

    @tb.task
    def write_results_to_csv(self, results_gap):
        rows = [
            {"name": name, "width": int(name.split("_")[1]), "gap": gap}
            for name, gap in results_gap.items()
        ]
        csv_path = Path("results.csv")
        with open(csv_path, mode="w", newline="") as csvfile:
            writer = csv.DictWriter(csvfile,
                                    fieldnames=["name", "width", "gap"])
            writer.writeheader()
            writer.writerows(rows)

        return csv_path


def workflow(runner):
    runner.run_workflow(Workflow(max_width=20))
