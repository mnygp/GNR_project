import matplotlib.pyplot as plt

import scienceplots
plt.style.use('science')

plt.rcParams['axes.grid'] = True

folder = "plots/"

k = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
e_k = [-10.00316365, -10.07809841, -10.06853679, -10.07162967,
       -10.06976777, -10.06998079, -10.0698612,  -10.06998801,
       -10.07030396, -10.07007983]

cut = [200, 240, 280, 320, 360, 400, 440, 480, 520, 560, 600]
e_cut = [-6.62525427, -8.35617688, -9.38751883, -9.88214847,
         -10.11034076, -10.22090821, -10.27173915, -10.29501596, -10.30432329,
         -10.30737315, -10.30854873]

vac = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14]
e_vac = [-9.80283683, -10.05942359, -10.0708293, -10.06998801,
         -10.06936937, -10.06863834, -10.07102616, -10.07032949, -10.06982983,
         -10.07052531, -10.06920361, -10.06986037, -10.06869701]


fig_k, ax_k = plt.subplots()
fig_k.set_size_inches(6, 4)
ax_k.plot(range(1, 11), e_k, 'o-')
ax_k.set_xlabel('k-points')
ax_k.set_ylabel('Energy (eV)')
ax_k.set_title('Convergence with k-points')
fig_k.savefig('../' + folder + 'AGNR_S_3_1_k_conv.png', dpi=500)

fig_cut, ax_cut = plt.subplots()
fig_cut.set_size_inches(6, 4)
ax_cut.plot(range(200, 601, 40), e_cut, 'o-')
ax_cut.set_xlabel('Plane-wave cutoff (eV)')
ax_cut.set_ylabel('Energy (eV)')
ax_cut.set_title('Convergence with plane-wave cutoff')
fig_cut.savefig('../' + folder + 'AGNR_S_3_1_cut_conv.png', dpi=500)

fig_vac, ax_vac = plt.subplots()
fig_vac.set_size_inches(6, 4)
ax_vac.plot(range(2, 15), e_vac, 'o-')
ax_vac.set_xlabel('Vacuum (Å)')
ax_vac.set_ylabel('Energy (eV)')
ax_vac.set_title('Convergence with vacuum')
fig_vac.savefig('../' + folder + 'AGNR_S_3_1_vac_conv.png', dpi=500)
