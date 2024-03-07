# Comparing the 4th pegs derivaties calculated by

# 1) FD derivatie before any other processing/bandpassing
# 2) Spectral derivative before any other processing
# 3) FD derivative of cut signal
# 4) Spectral derivative of cut signal

import numpy as np
import matplotlib.pyplot as plt
import obspy
from wetools import obspy_gen_mpl
stn = 'MDJ'
fig, ax = plt.subplots(figsize=(12,7))

ax.axvline(198.81, color='k', linestyle=':')

# Load signals differentiated before bandpassing:


# ------------  Local FD derivative:   ------------
fd_bef_path = "FD_local_deriv"
fd_bef_zacc_tr = obspy.read(f"{fd_bef_path}/Zacc/{stn}_proc.sac")[0]
fd_bef_grav_tr = obspy.read(f"{fd_bef_path}/Grav/{stn}_proc.sac")[0]

t_df_b, fd_bef_zacc = obspy_gen_mpl(fd_bef_zacc_tr)
t_df_b, fd_bef_grav = obspy_gen_mpl(fd_bef_grav_tr)
#ax.plot(t_df_b, fd_bef_zacc+fd_bef_grav, 'k')



# ------------  Local FD derivative NOT CUT:   ------------
fd_bef_path = "FD_local_deriv"
fd_bef_zacc_tr_notcut = obspy.read(f"{fd_bef_path}/Zacc/{stn}_proc_notcut.sac")[0]
fd_bef_grav_tr_notcut = obspy.read(f"{fd_bef_path}/Grav/{stn}_proc_notcut.sac")[0]

t_df_b_nc, fd_bef_zacc_notcut = obspy_gen_mpl(fd_bef_zacc_tr_notcut)
t_df_b_nc, fd_bef_grav_notcut = obspy_gen_mpl(fd_bef_grav_tr_notcut)
ax.plot(t_df_b_nc, fd_bef_zacc_notcut+fd_bef_grav_notcut, 'b')




# ------------  Spectral derivative:   ------------
spec_bef_path = "spectral_deriv"
spec_bef_zacc_tr = obspy.read(f"{spec_bef_path}/Zacc/{stn}_proc.sac")[0]
spec_bef_grav_tr = obspy.read(f"{spec_bef_path}/Grav/{stn}_proc.sac")[0]

t_spec_b, spec_bef_zacc = obspy_gen_mpl(spec_bef_zacc_tr)
t_spec_b, spec_bef_grav = obspy_gen_mpl(spec_bef_grav_tr)
#ax.plot(t_spec_b, spec_bef_zacc+spec_bef_grav, 'g')




# ------------  Spectral derivative NOT CUT:   ------------
spec_bef_path = "spectral_deriv"
spec_bef_zacc_tr_notcut = obspy.read(f"{spec_bef_path}/Zacc/{stn}_proc_notcut.sac")[0]
spec_bef_grav_tr_notcut = obspy.read(f"{spec_bef_path}/Grav/{stn}_proc_notcut.sac")[0]

t_spec_b_nc, spec_bef_zacc_notcut = obspy_gen_mpl(spec_bef_zacc_tr_notcut)
t_spec_b_nc, spec_bef_grav_notcut = obspy_gen_mpl(spec_bef_grav_tr_notcut)
ax.plot(t_spec_b_nc, spec_bef_zacc_notcut+spec_bef_grav_notcut, 'r')




# ------------ Spline derivative NOT CUT:   ------------
for n in [3]:
    spline_bef_path = f"spline_deriv/n{n}"
    spline_bef_zacc_tr_notcut = obspy.read(f"{spline_bef_path}/Zacc/{stn}_proc_notcut.sac")[0]
    spline_bef_grav_tr_notcut = obspy.read(f"{spline_bef_path}/Grav/{stn}_proc_notcut.sac")[0]

    t_spline_b_nc, spline_bef_zacc_notcut = obspy_gen_mpl(spline_bef_zacc_tr_notcut)
    t_spline_b_nc, spline_bef_grav_notcut = obspy_gen_mpl(spline_bef_grav_tr_notcut)
    ax.plot(t_spline_b_nc, spline_bef_zacc_notcut+spline_bef_grav_notcut)







ax.axvline(0.85*t_spec_b_nc[-1], alpha=0.15)
ax.axvline(0.15*t_spec_b_nc[-1], alpha=0.15)

#pm = 1e-11
#ax.set_xlim([0, 220])
#ax.set_ylim([-pm, pm])

ax.legend(['P arrival', 'FD', 'Spectral', 'Spline (n=3)', 'Tukey taper corners'])

plt.show()