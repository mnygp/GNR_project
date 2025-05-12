from functions.generate_ribbons import edge_state_ribbon
from functions.bandstructure import LDOS_from_file, LDOS

import matplotlib.pyplot as plt
import numpy as np

from pathlib import Path

from ase.io import read, write
from ase import Atoms

from gpaw import GPAW, PW, FermiDirac, restart

import scienceplots  # noqa: F401
plt.style.use('science')


def generate_ribbon(repeat, n, m) -> Path:
    length = 2*n + 2*(m - 1)
    # The edge state tag is 10
    ribbon = edge_state_ribbon(7, 'S', n, m, clean_length=length*repeat,
                               repeat=repeat, vac=5, tag=True)
    write(f'7-AGNR-repeat-{repeat}.xyz', ribbon)
    return Path(f'7-AGNR-repeat-{repeat}.xyz')


def tag_atoms(atoms_path: Path, m: int) -> Path:
    atoms = read(atoms_path.resolve())
    width = atoms.get_cell()[0, 0]
    height = atoms.get_cell()[1, 1]
    length = atoms.get_cell()[2, 2]

    middle_clean = (width/2, height/2, 3*length/4)
    middle_edge = (width/2, height/2, length/4)

    positions = atoms.positions

    # Closest atoms to the middle of the clean part
    clean_edge = np.array([np.linalg.norm(pos - middle_clean)
                           for pos in positions])
    clean_edge_index = np.argmin(clean_edge)
    atoms[clean_edge_index].tag = 20

    # Closest atoms to the middle of the serated part
    edge_state = np.array([np.linalg.norm(pos - middle_edge)
                           for pos in positions])
    edge_state_index = np.argmin(edge_state)
    atoms[edge_state_index].tag = 30

    write(f'7-AGNR-m-{m}-tagged.xyz', atoms)
    return Path(f'7-AGNR-m-{m}-tagged.xyz')


def gs_calculation(atoms_path: Path, params: dict,
                   file_name: str) -> Path:

    atoms = read(atoms_path)
    if params.get('PW_cut') is not None:
        atoms.calc = GPAW(mode=PW(params.get('PW_cut')),
                          xc=params['func'],
                          kpts={'size': (1, 1, 1)},
                          occupations=FermiDirac(0.01),
                          txt=file_name)
    else:
        atoms.calc = GPAW(mode='lcao',
                          basis=params['basis'],
                          xc=params['func'],
                          kpts={'size': (1, 1, 1)},
                          occupations=FermiDirac(0.01),
                          txt=file_name)

    atoms.get_potential_energy()
    atoms.calc.write(file_name, mode='all')
    return Path(file_name)


def LDOS_func(gpaw_path: Path, tag_number: int) -> dict:
    atoms, calc = restart(gpaw_path)
    ef = calc.get_fermi_level()
    print(f'Fermi energy {ef}')
    tags = atoms.get_tags()
    site_index = np.where(tags == tag_number)[0]
    energies, ldos = LDOS_from_file(gpaw_path, site_index, npoints=1200)
    return {'energies': energies.tolist(), 'ldos': ldos.tolist()}


def gs_and_LDOS(atoms_path: Path, tag_number: int):
    atoms = read(atoms_path)
    ldos_params = {'basis': 'szp(dzp)', 'func': 'PBE'}
    energies, ldos = LDOS(atoms, ldos_params, tag_number, k=1, npoints=1200)
    return {'energies': energies.tolist(), 'ldos': ldos.tolist()}


def insert_OH(atoms_path: Path) -> Path:
    OH_bond = 0.96
    CO_bond = 1.43
    atoms = read(atoms_path)
    symbols = atoms.get_chemical_symbols()
    site_index = np.where(atoms.get_tags() == 10)[0]

    # Add the oxygen atom
    H_atom_indices = [i for i, s in enumerate(symbols) if s == 'H']
    H_atom_distances = [np.linalg.norm(atoms.get_positions()[i]
                                       - atoms.get_positions()[site_index])
                        for i in H_atom_indices]
    closest_H_index = H_atom_indices[np.argmin(H_atom_distances)]
    atoms[closest_H_index].symbol = 'O'
    atoms[closest_H_index].tag = 40
    O_index = closest_H_index

    # Move it to the correct distance
    direction = atoms.positions[O_index] - atoms.positions[site_index]
    new_direction = direction/np.linalg.norm(direction)*CO_bond
    atoms.positions[O_index] = atoms.positions[site_index] + new_direction

    O_pos = atoms.get_positions()[O_index]
    # The *0.95 is to compensate that it get moved up a bit
    new_H_pos = np.array([O_pos[0] + np.sin(np.pi/6)*OH_bond,
                          O_pos[1] + 0.3,
                          O_pos[2] + np.cos(np.pi/6)*OH_bond])

    atoms += Atoms('H', positions=[new_H_pos])

    write('7-AGNR-OH.xyz', atoms)
    return Path('7-AGNR-OH.xyz')


def plot_ldos(energies, ldos, title: str, filename: str):
    plt.figure()
    plt.plot(energies, ldos)
    plt.title(title)
    plt.xlabel('Energy [eV]')
    plt.xlim(-3, 3)
    plt.ylim(0, 1.8)
    plt.ylabel('LDOS')
    plt.tight_layout()
    plt.savefig(filename, dpi=500)
    plt.close()
