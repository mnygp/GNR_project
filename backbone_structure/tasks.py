import taskblaster as tb

from ase.build import graphene_nanoribbon
from ase.io import write, read

from pathlib import Path
import csv

from functions.relax import relax
from functions.bandstructure import get_gap


# This is a taskblaster workflow that generates a series of workflows for
# generating armchair graphene nanoribbons (AGNR) of different widths.
@tb.dynamical_workflow_generator_task
def generate_wfs_from_list(input):
    for i in range(2, input+1):
        wf = SubWorkflow(width=i)
        name = f'AGNR_{i}'
        yield name, wf


# This is a taskblaster workflow that generates a series of workflows for
# generating armchair graphene nanoribbons (AGNR) of different widths.
@tb.workflow
class SubWorkflow:
    width = tb.var()

    @tb.task
    def generate_ribbon_task(self):
        return tb.node('generate_ribbon', width=self.width)

    @tb.task
    def gap_pre_relaxation_task(self):
        return tb.node('calculate_gap', in_file=self.generate_ribbon_task)

    @tb.task
    def relax_ribbon_task(self):
        return tb.node('relax_ribbon', in_file=self.generate_ribbon_task,
                       out_file=f'relaxed_ribbon_{self.width}.xyz')

    @tb.task
    def gap_post_relaxation_task(self):
        return tb.node('calculate_gap', in_file=self.relax_ribbon_task)

    @tb.task
    def return_dict_task(self):
        return tb.node('return_dict',
                       width=self.width,
                       pre_gap=self.gap_pre_relaxation_task,
                       post_gap=self.gap_post_relaxation_task)


# The functions called in the tasks above are defined here.
def generate_ribbon(width) -> Path:
    ribbon = graphene_nanoribbon(n=width/2.0, m=2, type='armchair',
                                 vacuum=6.0, saturated=True)
    write(f'AGNR_{width}.xyz', ribbon)
    return Path(f'AGNR_{width}.xyz')


def relax_ribbon(in_file: Path, out_file: str) -> Path:
    atoms = read(in_file)
    relaxed_ribbon = relax(atoms, 'relaxed_ribbon.txt', PW_toggle=True,
                           params={'func': 'PBE', 'PW_cut': 600}, k=4)
    write(out_file, relaxed_ribbon)
    return Path(out_file)


def calculate_gap(in_file: Path) -> float:
    ribbon = read(in_file)
    gap = get_gap(ribbon, params={'func': 'PBE', 'PW_cut': 600}, k=6,
                  filename='gap.txt')
    return gap


def return_dict(width, pre_gap, post_gap) -> dict:
    return {'width': width, 'pre_gap': pre_gap, 'post_gap': post_gap}


def write_results_to_csv(results) -> Path:
    rows = []
    for name, d in results.items():
        width = d['width']
        pre_gap = d['pre_gap']
        post_gap = d['post_gap']

        rows.append({
            "name": name,
            "width": width,
            "pre_relax_gap": pre_gap,
            "post_relax_gap": post_gap,
        })

    csv_path = Path("results.csv")
    with open(csv_path, mode="w", newline="") as csvfile:
        writer = csv.DictWriter(csvfile,
                                fieldnames=["name",
                                            "width",
                                            "pre_relax_gap",
                                            "post_relax_gap"])
        writer.writeheader()
        writer.writerows(rows)

    return csv_path
