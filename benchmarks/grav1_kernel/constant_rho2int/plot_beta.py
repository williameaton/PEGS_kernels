# The idea here will be to try and calculate the convolution of the gravity kernels directly
import ensightreader
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from load_proc_data import ProcData
from matplotlib.gridspec import GridSpec

def norm(x):
    return x/np.max(np.abs(x))

G = 6.6723e-11
X = []; Y= []; Z = []; V=[]


fpath = '/Users/eaton/REGS1/'
nprocs = 96
for i_proc in range(nprocs):
    print(i_proc)
    DProc = ProcData(proc=i_proc, fpath=fpath)


    DProc.load_var('gravbeta11')
    g1 = DProc.data["gravbeta11"]

    mask = np.logical_and(DProc.coord[:,2] < 0.000000000005, DProc.coord[:,2] > -0.000000000005)
    x = DProc.coord[mask,0]
    y = DProc.coord[mask,1]
    z = DProc.coord[mask,2]
    v = g1[mask]
    X = np.append(X,x[:])
    Y = np.append(Y,y[:])
    Z = np.append(Z,z[:])
    V = np.append(V,v[:])


fig, ax = plt.subplots()


R = 6371000
Rcmb = (1-2891/6371)*R


rplot = np.linspace(Rcmb, R, 10000)
beta  = (2/3)*np.pi*G * (3*(R**2) - rplot**2  - 2*(Rcmb**3)/rplot )


ax.plot(rplot, beta )


r = R * (X**2 + Y**2 + Z**2)**0.5
ax.plot(r, V/30, 'x')

plt.show()