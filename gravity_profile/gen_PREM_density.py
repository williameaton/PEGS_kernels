import matplotlib.pyplot as plt
import numpy as np

G = 6.67e-11

def get_PREM_density(r, R=6371000):

    x = r/R

    # Inner core:
    icmask = np.logical_and(r>=0, r<=1221500)
    x_ic =x[icmask]
    r_ic =r[icmask]
    rho_ic = 13.0885 - 8.8381* (x_ic**2)

    # Outer core:
    ocmask = np.logical_and(r>=1221500, r<=3480000)
    x_oc =x[ocmask]
    r_oc =r[ocmask]
    rho_oc = 12.5815 - 1.2638*x_oc - 3.6426*(x_oc**2) - 5.5281*(x_oc**3)

    # Lower mantle
    lmmask = np.logical_and(r>=3480000, r<=5701000)
    x_lm =x[lmmask]
    r_lm =r[lmmask]
    rho_lm = 7.9565 - 6.4761*x_lm + 5.5283*(x_lm**2) - 3.0807*(x_lm**3)

    # Transition zone 1
    tz1mask = np.logical_and(r>=5701000, r<=5771000)
    x_tz1 =x[tz1mask]
    r_tz1 =r[tz1mask]
    rho_tz1 = 5.3197 - 1.4836*x_tz1

    # Transition zone 2
    tz2mask = np.logical_and(r>=5771000, r<=5971000)
    x_tz2 =x[tz2mask]
    r_tz2 =r[tz2mask]
    rho_tz2 = 11.2494 -8.0298*x_tz2

    # Transition zone 3
    tz3mask = np.logical_and(r>=5971000, r<=6151000)
    x_tz3 =x[tz3mask]
    r_tz3 =r[tz3mask]
    rho_tz3 = 7.1089 -3.8045*x_tz3

    # Low velocity zone/ LID
    lvzmask = np.logical_and(r>=6151000, r<=6346600)
    x_lvz =x[lvzmask]
    r_lvz =r[lvzmask]
    rho_lvz = 2.6910 + 0.6924*x_lvz

    # Lower crust
    lcmask = np.logical_and(r >= 6346600, r <= 6356000)
    x_lc = x[lcmask]
    r_lc = r[lcmask]
    rho_lc = 2.9000 + x_lc*0

    # Upper crust
    ucmask = np.logical_and(r >= 6356000, r <= 6371000)
    x_uc = x[ucmask]
    r_uc = r[ucmask]
    rho_uc = 2.6000 + x_uc*0


    radius = np.concatenate((r_ic, r_oc, r_lm, r_tz1, r_tz2,  r_tz3, r_lvz, r_lc, r_uc))
    density = np.concatenate((rho_ic, rho_oc, rho_lm, rho_tz1, rho_tz2, rho_tz3, rho_lvz, rho_lc, rho_uc))

    density*=1000
    return radius, density




rr = np.linspace(0, 6371000, 50000)

dr = np.mean(rr[1:]-rr[:-1])


mp = rr[:-1] + dr/2
mp, density = get_PREM_density(mp)

#density = density*0 + 2600

phi = mp*0
g   = mp*0

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



fig, ax = plt.subplots(3)
ax[0].plot(mp, density)
ax[1].plot(mp, phi)



ax[1].plot(mp, phi)
ax[2].plot(mp, g)




plt.show()