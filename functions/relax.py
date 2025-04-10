from ase.optimize.bfgs import BFGS
from ase.filters import UnitCellFilter
from ase.io import Trajectory
from ase import Atoms

from gpaw import GPAW, PW

from typing import TypedDict, NotRequired

"""
This module contains functions to relax a structure using GPAW.
It provides functions to perform a single relaxation step or multiple
relaxation steps with different k-point meshes or basis sets.
The `relax` function performs a single relaxation step, while the
`multi_step_relax` function allows for multiple relaxation steps
with different parameters.
The `trajectory_file_name` function generates a filename for the
trajectory based on the provided parameters.

The parameters dictionary for the relaxation can include:
- `func`: The exchange-correlation functional to use (e.g., 'PBE').
- `PW_cut`: The plane-wave cutoff energy (in eV) for the GPAW calculator.
- `basis`: The basis set to use (e.g., 'dzp')
- `basis_list`: A list of basis sets to use for multiple relaxation steps.
"""


class RelaxParams(TypedDict):
    func: str
    PW_cut: NotRequired[float]
    basis: NotRequired[str]
    basis_list: NotRequired[list[str]]


def relax(structure: Atoms, filename: str, PW_toggle: bool,
          params: RelaxParams, k: int = 6,
          traj_name: str | None = None) -> Atoms:

    func = params['func']
    if PW_toggle:
        PW_cut = params.get('PW_cut')
        assert PW_cut is not None, "PW_cut is required when PW_toggle is True"
        calc = GPAW(mode=PW(PW_cut),
                    xc=func,
                    kpts={'size': (1, 1, k)},
                    txt=filename)
    else:
        basis = params.get('basis')
        assert basis is not None, "basis is required when PW_toggle is False"
        calc = GPAW(mode='lcao',
                    basis=basis,
                    xc=func,
                    kpts={'size': (1, 1, k)},
                    txt=filename)

    structure.calc = calc
    uf = UnitCellFilter(structure)

    if traj_name is not None:
        traject = Trajectory(traj_name+'.traj', 'w', uf)
        relax = BFGS(uf, trajectory=traject)
    else:
        relax = BFGS(uf)
    relax.run(fmax=0.01)

    return structure


def multi_step_relax(structure: Atoms, filename: str, PW_toggle: bool,
                     params: RelaxParams, k_arr: list[int],
                     traj_name: str | None = None) -> Atoms:

    if type(params['basis']) is list:
        if len(params['basis']) > 1 and len(params['basis']) != len(k_arr):
            raise ValueError("If multiple basis sets are provided," +
                             " the number of basis sets must match" +
                             " the number of k-points.")
        basis_list = params['basis']

    for i, k in enumerate(k_arr):
        if PW_toggle:
            traj_file_name = trajectory_file_name(traj_name, k,
                                                  PW_cut=params['PW_cut'])
        else:
            traj_file_name = trajectory_file_name(traj_name, k,
                                                  basis=params['basis'])

        if type(params['basis']) is list:
            params['basis'] = basis_list[i]

        structure = relax(structure, filename, PW_toggle, params, k,
                          traj_name=traj_file_name)

    return structure


def trajectory_file_name(traj_name: str | None, k: int,
                         basis: str | None = None, PW_cut: float | None = None
                         ) -> str | None:

    if traj_name is not None:
        if basis is None and PW_cut is not None:
            traj_file_name = traj_name + f'_k_{k}_{PW_cut}'
        elif basis is not None and PW_cut is None:
            traj_file_name = traj_name + f'_k_{k}_{basis}'
    else:
        traj_file_name = None

    return traj_file_name
