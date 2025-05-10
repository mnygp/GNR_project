import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d
from scipy.interpolate import interp1d

import scienceplots  # noqa: F401
plt.style.use('science')


energy_common = np.linspace(-2, 2, 400)
m_vals = [1, 2, 3, 4, 5]


for n in [1, 3]:
    ldos_all = []
    for m in range(1, 6):
        data = np.genfromtxt(f'tree/wfs_n_{n}/m_{m}/' +
                             f'write_csv_task/ldos_n_{n}_m_{m}.csv',
                             dtype=float, delimiter=',', skip_header=1)
        energies = data[:, 0]
        ldos = data[:, 1]

        interp = interp1d(energies, ldos, bounds_error=False, fill_value=0.0)
        ldos_interp = interp(energy_common)

        ldos_log = np.log(ldos_interp + 1)

        ldos_smooth = gaussian_filter1d(ldos_log, sigma=0.01)

        ldos_all.append(ldos_smooth)

    ldos_all = np.column_stack(ldos_all)

    plt.figure(figsize=(5, 4))
    plt.pcolormesh(m_vals, energy_common, ldos_all,
                   shading='auto', cmap='Reds')
    plt.colorbar(label='log(LDOS+1)')
    plt.xlabel('$m$')
    plt.ylabel('Energy (eV)')
    plt.text(4.7, 0.55, f'$n = {n}$', size='large')
    if (n == 1):
        plt.text(-0.2, 1.8, '\\textbf{a)}', size='x-large')
    else:
        plt.text(-0.2, 1.8, '\\textbf{b)}', size='x-large')
    plt.tight_layout()
    plt.savefig(f'gap_plot_n_{n}', dpi=500)
