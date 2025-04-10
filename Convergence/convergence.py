import numpy as np
from gpaw import GPAW, PW
from ase.parallel import parprint
import functions.generate_ribbons as gr

convergence_criterion = 0.01

folder = "convergence_files/"

cut_arr = np.array([])
energy_cut = np.array([])
k_arr = np.array([])
energy_k = np.array([])
vac_arr = np.array([])
energy_vac = np.array([])


cut = 200
diff = 200
while diff > convergence_criterion:

    ribbon = gr.generate_ribbon(3, 'S', 1, 1, vac=8)
    n_atoms = len(ribbon.positions)
    ribbon_saturated = gr.saturate_edges(ribbon)

    calc = GPAW(mode=PW(cut),
                xc='PBE',
                kpts={'size': (1, 1, 8)},
                txt=folder + f'AGNR_3_S_1_1_cut{cut}.txt')

    ribbon_saturated.calc = calc
    energy = ribbon_saturated.get_potential_energy()

    cut_arr = np.append(cut_arr, cut)
    energy_cut = np.append(energy_cut, energy/n_atoms)

    if len(energy_cut) > 1:
        diff = abs(energy_cut[-1] - energy_cut[-2])

    cut += 40
    parprint(f"Plane-wave cutoff {cut} with energy={energy/n_atoms:.3f}" +
             f" and diff={diff:.3f}")

if np.isnan(cut_arr[-1]) or cut_arr[-1] <= 0:
    raise ValueError(f"Invalid energy_cut value: {cut_arr[-1]}")

k = 1
diff = 200
while diff > convergence_criterion:

    ribbon = gr.generate_ribbon(3, 'S', 1, 1, vac=8)
    n_atoms = len(ribbon.positions)
    ribbon_saturated = gr.saturate_edges(ribbon)

    calc = GPAW(mode=PW(cut_arr[-1]),
                xc='PBE',
                kpts={'size': (1, 1, 8)},
                txt=folder + f'AGNR_3_S_1_1_k{k}.txt')

    ribbon_saturated.calc = calc
    energy = ribbon_saturated.get_potential_energy()

    k_arr = np.append(k_arr, k)
    energy_k = np.append(energy_k, energy/n_atoms)

    if len(energy_k) > 1:
        diff = abs(energy_k[-1] - energy_k[-2])

    k += 1
    parprint(f"k-points {k} with energy={energy/n_atoms:.3f}" +
             f" and diff={diff:.3f}")

if np.isnan(k_arr[-1]) or (k_arr[-1] % 1) != 0:
    raise ValueError(f"Invalid k-value value: {k_arr[-1]}")

vac = 3
diff = 200
while diff > convergence_criterion:

    ribbon = gr.generate_ribbon(3, 'S', 1, 1, vac=vac)
    n_atoms = len(ribbon.positions)
    ribbon_saturated = gr.saturate_edges(ribbon)

    calc = GPAW(mode=PW(cut_arr[-1]),
                xc='PBE',
                kpts={'size': (1, 1, int(k_arr[-1]))},
                txt=folder + f'AGNR_3_S_1_1_vac{vac}.txt')

    ribbon_saturated.calc = calc
    energy = ribbon.get_potential_energy()

    vac_arr = np.append(vac_arr, vac)
    energy_vac = np.append(energy_vac, energy/n_atoms)

    if len(energy_vac) > 1:
        diff = abs(energy_vac[-1] - energy_vac[-2])

    vac += 1
    parprint(f"Vacuum {vac}Å with energy={energy/n_atoms:.3f}" +
             f" and diff={diff:.3f}")


parprint("Convergence tests:")
parprint("k-points")
parprint(k_arr)
parprint(energy_k)
parprint("Plane-wave cutoff")
parprint(cut_arr)
parprint(energy_cut)
parprint("Vacuum 2Å to 14Å")
parprint(vac_arr)
parprint(energy_vac)
