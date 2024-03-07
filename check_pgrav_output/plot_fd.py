import numpy as np
import matplotlib.pyplot as plt

stn = 'MDJ'
net = 'IC'
dir = '../Tohoku/forward_data_and_plots/elastic_test/OUTPUT_FILES_ELASTIC/'

fig, ax = plt.subplots(2)

G = []
for depth in [0, 2, 6, 8]:
    gtmp = np.loadtxt(f"{dir}/{stn}_{depth}.{net}.MXG.sem.ascii")
    G.append(gtmp[:,1])

    ax[0].plot(gtmp[:,0], gtmp[:,1])


G = np.array(G)
# FD:
grav=  (-G[0,:] + 8*G[1,:] - 8*G[2,:] + G[3,:])/(12*2000)

ax[1].plot(gtmp[:,0], grav)



# Load the outputted value:
gg = np.loadtxt(f"{dir}/{stn}_{depth}.{net}.MXZC.PGRAV.sem.ascii")
ax[1].plot(gg[:,0], gg[:,1])



plt.show()
