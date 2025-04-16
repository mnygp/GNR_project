from functions.generate_ribbons import generate_ribbon, saturate_edges
from functions.relax import relax, multi_step_relax
from ase.io import read, write
from pathlib import Path


def ribbon_path(N: int, identifier: str, n: int, m: float | int):
    ribbon = generate_ribbon(N, identifier, n, m)
    saturated_ribbon = saturate_edges(ribbon, "H")
    write(f'{N}-AGNR-{identifier}.xyz', saturated_ribbon)
    return Path(f'{N}-AGNR-{identifier}.xyz')


def single_PW_relax(atoms_path: Path, filename: str, params: dict):
    atoms = read(atoms_path)

    relaxed_ribbon = relax(atoms, filename, PW_toggle=True,
                           params=params, k=6, traj_name='relax_PW')
    write('relaxed_ribbon_PW.xyz', relaxed_ribbon)
    return Path('relaxed_ribbon_PW.xyz')


def single_LCAO_relax(atoms_path: Path, filename: str, params: dict):
    atoms = read(atoms_path)

    relaxed_ribbon = relax(atoms, filename, PW_toggle=False,
                           params=params, k=6, traj_name='relax_LCAO')
    write('relaxed_ribbon_LCAO.xyz', relaxed_ribbon)
    return Path('relaxed_ribbon_LCAO.xyz')


def multi_relaxation(atoms_path: Path, filename: str,
                     params: dict, k_arr: list[int]):
    atoms = read(atoms_path)

    relaxed_ribbon = multi_step_relax(atoms, filename, PW_toggle=False,
                                      params=params, k_arr=k_arr,
                                      traj_name='relax_LCAO_multi_step')
    write('relaxed_ribbon_LCAO_multi_step.xyz', relaxed_ribbon)
    return Path('relaxed_ribbon_LCAO_multi_step.xyz')
