from ase import Atoms
from gpaw import GPAW, PW, FermiDirac


def get_gap(ribbon: Atoms, kpts=4, kpts_path: int = 60,
            filename: str | None = None) -> float:

    # Set up GPAW calculator
    calc = GPAW(mode=PW(600),
                xc='PBE',
                kpts={'size': (1, 1, kpts)},
                occupations=FermiDirac(0.01),
                txt=filename,
                convergence={'bands': 'occupied'})

    # Attach calculator to the ribbon
    ribbon.calc = calc

    ribbon.get_potential_energy()

    homo, lumo = calc.get_homo_lumo()
    gap = lumo - homo

    return gap
