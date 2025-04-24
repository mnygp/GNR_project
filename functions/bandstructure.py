from ase import Atoms
from ase.io import write
from gpaw import GPAW, PW, FermiDirac
from functions.relax import RelaxParams


# For now it only uses PW mode to calculate the gap
def get_gap(ribbon: Atoms, params: RelaxParams, k: int = 6,
            filename: str | None = None) -> float:

    PW_cut = params.get('PW_cut')
    assert PW_cut is not None, "PW_cut is required"
    func = params['func']
    calc = GPAW(mode=PW(PW_cut),
                xc=func,
                kpts={'size': (1, 1, k)},
                occupations=FermiDirac(0.01),
                txt=filename)

    ribbon.calc = calc

    ribbon.get_potential_energy()

    homo, lumo = calc.get_homo_lumo()
    gap = lumo - homo

    return gap


def LDOS(ribbon: Atoms, calc_params: RelaxParams, site_index: int, k: int = 6):
    PW_cut = calc_params.get('PW_cut')
    assert PW_cut is not None, "PW_cut is required"
    func = calc_params['func']
    calc = GPAW(mode=PW(PW_cut),
                xc=func,
                kpts={'size': (1, 1, k)},
                occupations=FermiDirac(0.01),
                txt='output.txt')

    ribbon.calc = calc

    ribbon.get_potential_energy()
    write('LDOS_gs.gpw', mode='all')

    energies, ldos = calc.get_wigner_seitz_ldos(site_index)

    return energies, ldos
