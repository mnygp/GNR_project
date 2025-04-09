from ase.optimize.bfgs import BFGS
from ase.filters import UnitCellFilter
from gpaw import GPAW, PW
from ase import Atoms


def relax(structure: Atoms, filename: str, PW_toggle: bool,
          params: dict[str, str | float], k: int = 6) -> Atoms:
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
    relax = BFGS(uf)
    relax.run(fmax=0.01)

    return structure


def multi_step_relax(structure: Atoms, filename: str, PW_toggle: bool,
                     params: dict[str], k_arr: list[int]) -> Atoms:

    for k in k_arr:
        structure = relax(structure, filename+f"_k_{k}", PW_toggle, params, k)

    return structure
