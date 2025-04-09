from ase.optimize.bfgs import BFGS
from ase.filters import UnitCellFilter
from ase.io import Trajectory
from ase import Atoms

from gpaw import GPAW, PW


def relax(structure: Atoms, filename: str, PW_toggle: bool,
          params: dict[str, str | float], k: int = 6,
          traj_name: str | None = None) -> Atoms:
    func = params['func']
    if PW_toggle:
        PW_cut = params['PW_cut']
        calc = GPAW(mode=PW(PW_cut),
                    xc=func,
                    kpts=(1, 1, k),
                    txt=filename)
    else:
        basis = params['basis']
        calc = GPAW(mode='lcao',
                    basis=basis,
                    xc=func,
                    kpts=(1, 1, k),
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
                     params: dict[str], k_arr: list[int],
                     traj_name: str | None = None) -> Atoms:

    if traj_name is not None:
        for k in k_arr:
            traj_file_name = traj_name + f'_k_{k}'
            structure = relax(structure, traj_name, PW_toggle, params, k,
                              traj_name=traj_file_name)
    else:
        for k in k_arr:
            structure = relax(structure, filename+f'_k_{k}',
                              PW_toggle, params, k)

    return structure
