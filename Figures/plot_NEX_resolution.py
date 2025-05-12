# Plot surface resolution as a function of NEX number:
import matplotlib.pyplot as plt
import numpy as np
from wetools.plotting import setup_we_mpl
setup_we_mpl()

circ = 2 * np.pi * 6371   # km
ngll = 5
minspeed = 3  # km/s
nex = np.array([192, 256, 368, 384, 400, 512, 528, 576, 800])

length_res = circ/(nex * 4 * ngll)
time_res   = 2* length_res/minspeed

fig, ax = plt.subplots(1,2)

# 4 chunks go around Earth and 5 GLL points
ax[0].plot(nex, length_res)
ax[1].plot(nex, time_res)

ax[0].set_ylabel('Length scale [km]')
ax[1].set_ylabel('Time scale [s]')

ax[0].axvline(256)
ax[0].axvline(384)

plt.show()