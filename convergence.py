import numpy as np
from matplotlib import pyplot as plt
from gpaw import GPAW, PW
from ase.parallel import parprint
import sys
sys.path.append('.')
import functions.generate_ribbons as gr

folder = "convergence_files/"

# Create a 3-AGNR with 2 saturated edges on each side

energy_k = np.array([])
energy_cut = np.array([])
energy_vac = np.array([])

for k in range(1, 11):
    ribbon = gr.generate_ribbon(3, 'S', 1, 1)
    ribbon_saturated = gr.saturate_edges(ribbon)
    calc = GPAW(mode=PW(350),
                xc='PBE',
                kpts=(1, 1, k),
                txt=folder + f'AGNR_3_S_1_1_k{k}.txt')

    ribbon_saturated.calc = calc

    energy = ribbon_saturated.get_potential_energy()
    energy_k = np.append(energy_k, energy)
    parprint(f"k-points {k} done")


for cut in range(200, 601, 40):
    ribbon = gr.generate_ribbon(3, 'S', 1, 1)
    ribbon_saturated = gr.saturate_edges(ribbon)
    calc = GPAW(mode=PW(cut),
                xc='PBE',
                kpts=(1, 1, 8),
                txt=folder + f'AGNR_3_S_1_1_cut{cut}.txt')

    ribbon_saturated.calc = calc

    energy = ribbon_saturated.get_potential_energy()
    energy_cut = np.append(energy_cut, energy)
    parprint(f"Plane-wave cutoff {cut} done")

for vac in range(1, 15):
    ribbon = gr.generate_ribbon(3, 'S', 1, 1, vac=vac)
    ribbon_saturated = gr.saturate_edges(ribbon)
    calc = GPAW(mode=PW(350),
                xc='PBE',
                kpts=(1, 1, 8),
                txt=folder + f'AGNR_3_S_1_1_vac{vac}.txt')

    ribbon_saturated.calc = calc

    energy = ribbon.get_potential_energy()
    energy_vac = np.append(energy_vac, energy)
    parprint(f"Vacuum {vac}Å done")

parprint("Convergence tests:")
parprint("k-points 1 to 10")
parprint(energy_k)
parprint("Plane-wave cutoff 200 to 600 in steps of 40")
parprint(energy_cut)
parprint("Vacuum 1Å to 10Å")
parprint(energy_vac)
