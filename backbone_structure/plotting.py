import numpy as np
import matplotlib.pyplot as plt

data = np.genfromtxt('tree/write_csv_task/results.csv',
                     delimiter=',',
                     skip_header=1)

order = data[:, 1].argsort()

width = data[order, 1]
pre = data[order, 2]
post = data[order, 3]

plt.plot(width, pre, '-o', label="pre relaxation")
plt.plot(width, post, '-o', label="post relaxation")
plt.xlim(1,21)
plt.xticks(np.linspace(2, 20, 10, dtype=int))
plt.xlabel("Width")
plt.ylabel("Gap energy [eV]")
plt.tight_layout()
plt.legend()
plt.grid()
plt.savefig("backbone_structure.png", dpi=500)
