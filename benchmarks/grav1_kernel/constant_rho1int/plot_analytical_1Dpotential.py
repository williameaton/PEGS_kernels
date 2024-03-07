import numpy as np
import matplotlib.pyplot as plt
from copy import copy

R = 6371000
Rcmb = (1- (2891000/R))*R
pi = np.pi
G = 6.67e-11
k = 1

r = np.linspace(0, R, 10000)

pot = copy(r)*0

coremask = r < Rcmb
mantlemask = r > Rcmb


pot[coremask] = 2*pi*G*k *(R**2 - r[coremask]**2)
pot[mantlemask] = (2*pi*G*k/3) * (3*(R**2) - r[mantlemask]**2 - 2*(Rcmb**3)/r[mantlemask]   )


fig,ax = plt.subplots()


ax.plot(r, pot, 'k')

plt.show()