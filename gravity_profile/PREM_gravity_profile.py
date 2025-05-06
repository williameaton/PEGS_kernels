# Computes prem gravitational potential and gravitational acceleration by
# radial integration

import matplotlib.pyplot as plt
import numpy as np
from wetools.models import get_PREM_density
G = 6.67e-11

# Radius
rr = np.linspace(0, 6371000, 50000)

# Radial spacing
dr = np.mean(rr[1:]-rr[:-1])

# Midpoints
mp = rr[:-1] + dr/2

# Get density at midpoints
mp, density = get_PREM_density(mp)

# Initialise phi and g
phi = mp*0
g   = mp*0

# loop through midpoints
for i in range(len(mp)):

    r_below = mp[:i]
    rho_below = density[:i]
    r_above = mp[i:]
    rho_above = density[i:]

    # --------------- FOR PHI ---------------:
    # Stuff below:
    inte = np.sum(r_below**2 * dr * rho_below)
    below = inte/mp[i]

    # Stuff above:
    above = np.sum(r_above* dr * rho_above)
    phi[i] = - 4*np.pi*G *(above+below)

    # --------------- FOR G ---------------:
    # Stuff below:
    inte = np.sum((r_below ** 2)* dr * rho_below)
    g[i] = 4*np.pi*G * inte / (mp[i]**2)

# Plot
fig, ax = plt.subplots(3, sharex=True, figsize=(12,7))
ax[0].plot(mp, density)
ax[1].plot(mp, phi)
ax[2].plot(mp, g)

ax[0].set_title('PREM Density')
ax[1].set_title('PREM gravitational potential')
ax[2].set_title('PREM gravitational acceleration')
ax[-1].set_xlabel('Radius [km]')

plt.show()