import matplotlib.pyplot as plt

# print(plt.rcParams['font.sans-serif'])
# plt.rcParams['font.sans-serif'] = ['Computer Modern Sans Serif']

folder = "convergence_files/"

k = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
e_k = [-240.0759277,  -241.8743619,  -241.64488299, -241.71911204,
       -241.67442637, -241.67953907, -241.67666885, -241.67971234,
       -241.68729495, -241.68191597]

cut = [200, 240, 280, 320, 360, 400, 440, 480, 520, 560, 600]
e_cut = [-159.00610238, -200.54824502, -225.30045204, -237.17156331,
         -242.64817817, -245.30179695, -246.52173962, -247.08038294,
         -247.30375903, -247.37695552, -247.40516961]

vac = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14]
e_vac = [-162.73828827, -235.26808402, -241.42616613, -241.69990325,
         -241.67971234, -241.66486494, -241.64732009, -241.70462784,
         -241.68790778, -241.67591599, -241.69260739, -241.66088671,
         -241.67664887, -241.64872833]


fig_k, ax_k = plt.subplots()
ax_k.plot(range(1, 11), e_k, 'o-')
ax_k.set_xlabel('k-points')
ax_k.set_ylabel('Energy (eV)')
ax_k.set_title('Convergence with k-points')
ax_k.grid()
fig_k.tight_layout()
fig_k.savefig(folder + 'AGNR_S_3_1_k.png')

fig_cut, ax_cut = plt.subplots()
ax_cut.plot(range(200, 601, 40), e_cut, 'o-')
ax_cut.set_xlabel('Plane-wave cutoff (eV)')
ax_cut.set_ylabel('Energy (eV)')
ax_cut.set_title('Convergence with plane-wave cutoff')
ax_cut.grid()
fig_cut.tight_layout()
fig_cut.savefig(folder + 'AGNR_S_3_1_cut.png')

fig_vac, ax_vac = plt.subplots()
ax_vac.plot(range(1, 15), e_vac, 'o-')
ax_vac.set_xlabel('Vacuum (Å)')
ax_vac.set_ylabel('Energy (eV)')
ax_vac.set_title('Convergence with vacuum')
ax_vac.grid()
fig_vac.tight_layout()
fig_vac.savefig(folder + 'AGNR_S_3_1_vac.png')
