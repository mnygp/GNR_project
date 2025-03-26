from ase.optimize.bfgs import BFGS
from ase.filters import UnitCellFilter
from gpaw import GPAW, PW


def relax_PW(structure, filename, PW_cut=500, k_pts=6, func='PBE'):
    calc = GPAW(mode=PW(PW_cut),
                xc=func,
                kpts=(1, 1, k_pts),
                txt=filename)

    structure.set_calculator(calc)
    uf = UnitCellFilter(structure)
    relax = BFGS(uf)
    relax.run(fmax=0.01)

    return structure


def relax_LCAO(structure, filename, k_pts=6, func='PBE', basis='dzp'):
    calc = GPAW(mode='lcao',
                basis=basis,
                xc=func,
                kpts=(1, 1, k_pts),
                txt=filename)

    structure.set_calculator(calc)
    uf = UnitCellFilter(structure)
    relax = BFGS(uf)
    relax.run(fmax=0.01)

    return structure
