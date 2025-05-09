from functions.generate_ribbons import edge_state_ribbon
from functions.bandstructure import LDOS_from_file

import matplotlib.pyplot as plt
import numpy as np

from pathlib import Path

from ase.io import read, write

from gpaw import GPAW, FermiDirac, restart

import taskblaster as tb


@tb.dynamical_workflow_generator_task
def generate_wfs(m_max, n):
    for i in range(1, m_max+1):
        wf = SubWorkflow(n=n, m=i)
        name = f'm_{i}'
        yield name, wf


@tb.workflow
class SubWorkflow:
    n = tb.var()
    m = tb.var()

    @tb.task
    def generate_ribbon_task(self):
        return tb.node('generate_ribbon',
                       repeat=5,
                       n=self.n,
                       m=self.m)

    @tb.task
    def tag_atoms_task(self):
        return tb.node('tag_atoms',
                       atoms_path=self.generate_ribbon_task,
                       n=self.n,
                       m=self.m)

    @tb.task
    def gs_task(self):
        return tb.node('gs_calculation',
                       atoms_path=self.tag_atoms_task,
                       file_name=f'gs-n-{self.n}-m-{self.m}.gpw')

    @tb.task
    def ldos_edge_state_task(self):
        return tb.node('LDOS_func',
                       gpaw_path=self.gs_task,
                       tag_number=10)

    @tb.task
    def plot_ldos_edge_task(self):
        return tb.node('plot_ldos',
                       energies=self.ldos_edge_state_task['energies'],
                       ldos=self.ldos_edge_state_task['ldos'],
                       n=self.n,
                       m=self.m)


def generate_ribbon(repeat, n, m) -> Path:
    length = 2*n + 2*(m - 1)
    # The edge state tag is 10
    ribbon = edge_state_ribbon(7, 'S', n, m, clean_length=2*length,
                               repeat=repeat, vac=5, tag=True)
    write(f'7-AGNR-m-{m}.xyz', ribbon)
    return Path(f'7-AGNR-m-{m}.xyz')


def gs_calculation(atoms_path: Path, file_name: str) -> Path:

    atoms = read(atoms_path)
    calc = GPAW(mode='lcao',
                basis='szp(dzp)',
                xc='PBE',
                kpts={'size': (1, 1, 1), 'gamma': True},
                occupations=FermiDirac(0.01))

    atoms.calc = calc

    atoms.get_potential_energy()
    calc.get_fermi_level()
    atoms.calc.write(file_name, mode='all')
    return Path(file_name)


def tag_atoms(atoms_path: Path, n: int, m: int) -> Path:
    atoms = read(atoms_path.resolve())
    width = atoms.get_cell()[0, 0]
    height = atoms.get_cell()[1, 1]
    length = atoms.get_cell()[2, 2]

    # The whole cell can be divided into 7 parts 5 serrated and 2 clean
    middle_clean = (width/2, height/2, 6.0/7*length)
    middle_edge = (width/2, height/2, (5.0/7)/2*length)

    positions = atoms.positions

    # Closest atoms to the middle of the clean part
    clean_edge = np.array([np.linalg.norm(pos - middle_clean)
                           for pos in positions])
    clean_edge_index = np.argmin(clean_edge)
    atoms[clean_edge_index].tag = 20

    # Closest atoms to the middle of the edge state part
    edge_state = np.array([np.linalg.norm(pos - middle_edge)
                           for pos in positions])
    edge_state_index = np.argmin(edge_state)
    atoms[edge_state_index].tag = 30

    write(f'7-AGNR-m-{m}-tagged.xyz', atoms)
    return Path(f'7-AGNR-m-{m}-tagged.xyz')


def LDOS_func(gpaw_path: Path, tag_number: int) -> dict:
    atoms, calc = restart(gpaw_path)
    tags = atoms.get_tags()
    site_index = np.where(tags == tag_number)[0]
    energies, ldos = LDOS_from_file(gpaw_path, site_index, 1200)
    return {'energies': energies.tolist(), 'ldos': ldos.tolist()}


def plot_ldos(energies, ldos, n, m):
    plt.figure()
    plt.plot(energies, ldos)
    plt.title(f'LDOS for n={n} and m={m}')
    plt.xlabel('Energy [eV]')
    plt.xlim(-3, 3)
    plt.grid()
    plt.savefig(f'ldos_n_{n}_m_{m}.png', dpi=500)
    plt.close()
