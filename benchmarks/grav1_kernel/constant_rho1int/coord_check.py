import numpy as np
import matplotlib.pyplot as plt

from analytical import analytical
d = np.loadtxt('NEX96/slurm-10372219.out', usecols=[0])

print(np.shape(d))

dd = d

print(np.shape(dd))


ddmask = np.abs(dd[:,2]) < 0.001


fig, ax = plt.subplots()
ax.scatter(dd[ddmask,0], dd[ddmask,1],
           c=dd[ddmask,-1],
           cmap='RdBu_r', s=1)

plt.show()