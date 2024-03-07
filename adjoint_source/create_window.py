import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
import numpy as np

# Create some points between 0 and T:
def create_window_taper(T, nknots):

    time = np.zeros(nknots+2)
    ttmp = np.linspace(0.6*T, 0.9*T, nknots)

    time[1:-1] = ttmp
    time[-1] = T
    val = np.zeros(nknots+2)
    val[1:-1] = 1.0

    cs = CubicSpline(x=time, y=val, bc_type=((1, 0.0), (1, 0.0)) )

    return cs