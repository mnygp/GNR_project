from ase.parallel import parprint
from ase.io import write
from gpaw import GPAW, PW, FermiDirac

from ase.build import graphene_nanoribbon


for i in range(2, 6):
    ribbon = graphene_nanoribbon(n=(i+1)/2, m=2, type='armchair',
                                 vacuum=5.0, saturated=True)
    ribbon.pbc = (1, 1, 1)
    write(f'AGNR_{i}.xyz', ribbon)

    calc = GPAW(mode=PW(600),
                xc='PBE',
                kpts={'size': (1, 1, 6)},
                occupations=FermiDirac(0.01),
                txt='output.txt')
    ribbon.calc = calc
    ribbon.get_potential_energy()
    ef = calc.get_fermi_level()

    calc_bs = calc.fixed_density(kpts={'path': [(0, 0, 0), (0, 0, 0.5)],
                                       'npoints': 100},
                                 symmetry='off')

    bs = calc_bs.band_structure()
    bs.plot(filename=f'bandstructure_{i}.png', show=False, emin=-3, emax=3)
    parprint(f'{i+1} bandstructure done')
