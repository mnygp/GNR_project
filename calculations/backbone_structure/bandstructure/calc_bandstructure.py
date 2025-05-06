from ase.io import write
from ase.build import graphene_nanoribbon

from gpaw import GPAW, PW, FermiDirac

for i in range(2, 6):
    ribbon = graphene_nanoribbon(n=(i+1)/2, m=1, type='armchair',
                                 vacuum=5.0, saturated=True)
    ribbon.pbc = (1, 1, 1)
    write(f'AGNR_{i}.xyz', ribbon)

    calc = GPAW(mode=PW(600),
                xc='PBE',
                kpts={'size': (1, 1, 100)},
                occupations=FermiDirac(0.01),
                txt='output.txt')

    ribbon.calc = calc
    ribbon.get_potential_energy()
    ef = calc.get_fermi_level()

    bs = calc.band_structure()
    bs = bs.subtract_reference()
    bs.plot(filename=f'bandstructure_{i}.png', emin=-3, emax=3, show=False)
