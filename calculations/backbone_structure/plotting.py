import numpy as np
import matplotlib.pyplot as plt

import scienceplots  # noqa: F401
plt.style.use('science')

data = np.genfromtxt('tree/write_csv_task/results.csv',
                     delimiter=',',
                     skip_header=1)

order = data[:, 1].argsort()

width = data[order, 1]
pre = data[order, 2]
post = data[order, 3]

# N = 3p
w0, pre0, post0 = width[::3], pre[::3], post[::3]
# N = 3p + 1
w1, pre1, post1 = width[1::3], pre[1::3], post[1::3]
# N = 3p + 2
w2, pre2, post2 = width[2::3], pre[2::3], post[2::3]

fig, ax = plt.subplots(figsize=(5, 2.5))

ax.plot(width, post, '-o', color='black')
ax.plot(w0, post0, '-o', color='blue')
ax.plot(w1, post1, '-o', color='red')
ax.plot(w2, post2, '-o', color='green')
ax.set_xlim(1, 21)
ax.set_ylim(bottom=0)
ax.set_xticks(np.linspace(2, 20, 10, dtype=int))
ax.set_xlabel("Width")
ax.set_ylabel("Gap energy [eV]")
fig.tight_layout()
ax.grid()
fig.savefig("backbone_structure.png", dpi=500)
