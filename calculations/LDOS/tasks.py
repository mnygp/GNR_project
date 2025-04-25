from functions.generate_ribbons import edge_state_ribbon
from functions.relax import multi_step_relax
from functions.bandstructure import LDOS
from ase.io import read, write
from pathlib import Path
import taskblaster as tb
import numpy as np
import matplotlib.pyplot as plt


@tb.dynamical_workflow_generator_task
def generate_wfs_from_list(input):
    for i in range(1, input+1):
        wf = SubWorkflow(repeat=i)
        name = f'repeat_{i}'
        yield name, wf


@tb.workflow
class SubWorkflow:
    repeat = tb.var()

    @tb.task
    def generate_ribbon_task(self):
        return tb.node('generate_ribbon', repeat=self.repeat)

    @tb.task
    def ldos_pre_relaxation_task(self):
        ldos_params = {"func": "PBE", "basis": 'dzp'}
        return tb.node('LDOS_func', atoms_path=self.generate_ribbon_task,
                       params=ldos_params)

    @tb.task
    def plot_ldos_pre_relax_task(self):
        return tb.node('plot_ldos',
                       energies=self.ldos_pre_relaxation_task['energies'],
                       ldos=self.ldos_pre_relaxation_task['ldos'],
                       repeat=self.repeat)

    @tb.task
    def multi_relaxation_task(self):
        params = {"func": "PBE", "basis": ["szp(dzp)", "szp(dzp)", "dzp"]}
        k_arr = [1, 3, 6]
        return tb.node('multi_relaxation',
                       filename='calc_out_multi.txt',
                       atoms_path=self.generate_ribbon_task,
                       params=params,
                       k_arr=k_arr)

    @tb.task
    def ldos_post_relaxation_task(self):
        ldos_params = {"func": "PBE", "basis": 'dzp'}
        return tb.node('LDOS_func', atoms_path=self.multi_relaxation_task,
                       params=ldos_params)

    @tb.task
    def plot_ldos_post_relax_task(self):
        return tb.node('plot_ldos',
                       energies=self.ldos_post_relaxation_task['energies'],
                       ldos=self.ldos_post_relaxation_task['ldos'],
                       repeat=self.repeat)


def generate_ribbon(repeat) -> Path:
    n = 1
    m = 1
    length = 2*n + 2*(m - 1)
    ribbon = edge_state_ribbon(7, 'S', n, m, clean_length=length*repeat,
                               repeat=repeat, vac=5, tag=True)
    write(f'7-AGNR-repeat_{repeat}.xyz', ribbon)
    return Path(f'7-AGNR-repeat_{repeat}.xyz')


def multi_relaxation(atoms_path: Path, filename: str,
                     params: dict, k_arr: list[int]):
    atoms = read(atoms_path)

    relaxed_ribbon = multi_step_relax(atoms, filename, PW_toggle=False,
                                      params=params, k_arr=k_arr,
                                      traj_name='relax_LCAO_multi_step')
    write('relaxed_ribbon_LCAO_multi_step.xyz', relaxed_ribbon)
    return Path('relaxed_ribbon_LCAO_multi_step.xyz')


def LDOS_func(atoms_path: Path, params: dict):
    atoms = read(atoms_path)
    tags = atoms.get_tags()
    site_index = np.where(tags == 10)[0]
    energies, ldos = LDOS(atoms, params, site_index, npoints=800)
    return {'energies': energies.tolist(), 'ldos': ldos.tolist()}


def plot_ldos(energies, ldos, repeat):
    plt.figure()
    plt.plot(energies, ldos)
    plt.title(f'LDOS for {repeat} repetitions')
    plt.xlabel('Energy [eV]')
    plt.xlim(-2, 2)
    plt.grid()
    plt.savefig(f'ldos_repeat_{repeat}.png', dpi=500)
    plt.close()
