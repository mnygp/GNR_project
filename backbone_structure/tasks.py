import taskblaster as tb

from ase.build import graphene_nanoribbon
from ase.io import write, read

from pathlib import Path

from functions.relax import relax
from functions.bandstructure import get_gap


# This is a taskblaster workflow that generates a series of workflows for
# generating armchair graphene nanoribbons (AGNR) of different widths.
@tb.dynamical_workflow_generator_task
def generate_wfs_from_list(input):
    for i in range(4, input+1):
        wf = SubWorkflow(width=i)
        name = f'AGNR_{i}'
        yield name, wf


@tb.workflow
class SubWorkflow:
    width = tb.var()

    @tb.task
    def generate_ribbon_task(self):
        return tb.node('generate_ribbon', width=self.width)

    @tb.task
    def relax_ribbon_task(self):
        return tb.node('relax_ribbon', in_file=self.generate_ribbon_task,
                       out_file=f'relaxed_ribbon_{self.width}.xyz')

    @tb.task
    def calculate_gap_task(self):
        return tb.node('calculate_gap_task', in_file=self.relax_ribbon_task)


def generate_ribbon(width) -> Path:
    ribbon = graphene_nanoribbon(n=width/2.0, m=2, type='armchair',
                                 vacuum=8.0, saturated=True)
    write(f'AGNR_{width}.xyz', ribbon)
    return Path(f'AGNR_{width}.xyz')


def relax_ribbon(in_file: Path, out_file: str) -> Path:
    atoms = read(in_file)
    relaxed_ribbon = relax(atoms, 'relaxed_ribbon.txt', PW_toggle=True,
                           params={'func': 'PBE', 'PW_cut': 600}, k=4)
    write(out_file, relaxed_ribbon)
    return Path(out_file)


def calculate_gap_task(in_file: Path) -> float:
    ribbon = read(in_file)
    gap = get_gap(ribbon, params={'func': 'PBE', 'PW_cut': 600}, k=6,
                  filename='gap.txt')
    return gap
