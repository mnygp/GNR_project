from generate_ribbons import generate_ribbon, saturate_edges
from gpaw import GPAW, PW
import numpy as np

folder = "convergence_files/"

# Create a 3-AGNR with 2 saturated edges on each side
ribbon = generate_ribbon(3, 'S', 1, 1)
ribbon_saturated = saturate_edges(ribbon)

energy_k = np.array([])
energy_cut = np.array([])

for k in range(1, 11):
    calc = GPAW(mode=PW(350),
                xc='PBE',
                kpts=(1, 1, k),
                txt=folder + f'AGNR_S_3_1_k{k}.txt')

    ribbon_saturated.set_calculator(calc)

    energy = ribbon_saturated.get_potential_energy()
    energy_k = np.append(energy_k, energy)


for cut in range(300, 601, 50):
    calc = GPAW(mode=PW(cut),
                xc='PBE',
                kpts=(1, 1, 6),
                txt=folder + f'AGNR_S_3_1_cut{cut}.txt')

    ribbon_saturated.set_calculator(calc)

    energy = ribbon_saturated.get_potential_energy()
    energy_cut = np.append(energy_k, energy)