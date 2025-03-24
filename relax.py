from ase.optimize.bfgs import BFGS
from ase.filters import UnitCellFilter
from gpaw import GPAW, PW


def relax_PW(structure, PW_cut, k_pts, filename, func='PBE'):
    calc = GPAW(mode=PW(PW_cut),
                xc=func,
                kpts=(1, 1, k_pts),
                txt=filename)

    structure.set_calculator(calc)
    uf = UnitCellFilter(structure)
    relax = BFGS(uf)
    relax.run(fmax=0.05)

    return structure
