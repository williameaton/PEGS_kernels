import numpy as np
import matplotlib.pyplot as plt


analytical_color = 'cornflowerblue'
zhang_color = 'royalblue'


def reproduce_ZhangFig2(fpath):
    fig, ax = plt.subplots(3,2, figsize=(12,7))
    lw = 1
    for irow in range(3):
        stn = irow + 1
        # Analytical data # Hpp is 2nd column in file after time
        anal = np.loadtxt(f"{fpath}/Analytical/grstrain.RCVR{stn}.TOHOKU.dat")
        ax[irow,0].plot(anal[:,0], anal[:,1], color=analytical_color, linewidth=lw, label='Analytical Half Space')

        # Halfspace data QSSP:
        qssp = np.loadtxt(f"{fpath}/HalfSpace/RCVR{stn}_grav.h.dat")
        ax[irow,0].plot(qssp[:,0], qssp[:,1], color=zhang_color, linewidth=lw, label='QSSP (Zhang)')

        ax[irow,0].set_xlim([-10, 50])
        ax[irow,0].set_ylabel('Gravity Strain (x 1e-12)')

    for irow in range(3):
        stn = irow + 4
        # Analytical data # Hpp is 2nd column in file after time
        anal = np.loadtxt(f"{fpath}/Analytical/grstrain.RCVR{stn}.TOHOKU.dat")
        ax[irow,1].plot(anal[:,0], anal[:,1], color=analytical_color, linewidth=lw, label='Analytical Half Space')

        # Halfspace data QSSP:
        qssp = np.loadtxt(f"{fpath}/HalfSpace/RCVR{stn}_grav.h.dat")
        ax[irow,1].plot(qssp[:,0], qssp[:,1], color=zhang_color, linewidth=lw, label='QSSP (Zhang)')



        ax[irow,1].set_xlim([-20, 150])




    # Set Y lims:
    ax[0,0].set_ylim([-6, 1])
    ax[1,0].set_ylim([-8, 1])
    ax[2,0].set_ylim([-4, 1])

    ax[0,1].set_ylim([-30, 5])
    ax[1,1].set_ylim([-30,  5])
    ax[2,1].set_ylim([-35,  5])


    ax[2,0].set_xlabel('Time [s]')
    ax[2,1].set_xlabel('Time [s]')



    return fig, ax





def reproduce_ZhangFig3(fpath):

    fig, ax = plt.subplots(3,2, figsize=(12,7))
    lw = 1
    for irow in range(3):
        stn = irow + 1
        # Analytical data # Hpp is 2nd column in file after time
        anal = np.loadtxt(f"{fpath}/Analytical/grstrain.RCVR{stn}.TOHOKU.dat")
        ax[irow,0].plot(anal[:,0], anal[:,5], color=analytical_color, linewidth=lw, label='Analytical Half Space')

        # Halfspace data QSSP:
        qssp = np.loadtxt(f"{fpath}/HalfSpace/RCVR{stn}_grav.h.dat")
        ax[irow,0].plot(qssp[:,0], qssp[:,5], color=zhang_color, linewidth=lw, label='QSSP (Zhang)')

        ax[irow,0].set_xlim([-10, 50])
        ax[irow,0].set_ylabel('Gravity Strain (x 1e-12)')

    for irow in range(3):
        stn = irow + 4
        # Analytical data # Hpp is 2nd column in file after time
        anal = np.loadtxt(f"{fpath}/Analytical/grstrain.RCVR{stn}.TOHOKU.dat")
        ax[irow,1].plot(anal[:,0], anal[:,5], color=analytical_color, linewidth=lw, label='Analytical Half Space')

        # Halfspace data QSSP:
        qssp = np.loadtxt(f"{fpath}/HalfSpace/RCVR{stn}_grav.h.dat")
        ax[irow,1].plot(qssp[:,0], qssp[:,5], color=zhang_color, linewidth=lw, label='QSSP (Zhang)')



        ax[irow,1].set_xlim([-20, 150])




    # Set Y lims:
    ax[0,0].set_ylim([-1, 4])
    ax[1,0].set_ylim([-2, 6])
    ax[2,0].set_ylim([-1, 3])

    ax[0,1].set_ylim([-4, 6])
    ax[1,1].set_ylim([-5, 5])
    ax[2,1].set_ylim([-4, 6])


    ax[2,0].set_xlabel('Time [s]')
    ax[2,1].set_xlabel('Time [s]')

    return fig, ax






def plot_zhang(fpath, index):

    fig, ax = plt.subplots(3,2, figsize=(12,7))
    lw = 1
    for irow in range(3):
        stn = irow + 1
        anal = np.loadtxt(f"{fpath}/Analytical/grstrain.RCVR{stn}.TOHOKU.dat")
        ax[irow,0].plot(anal[:,0], anal[:,index], color=analytical_color, linewidth=lw, label='Analytical Half Space')

        # Halfspace data QSSP:
        qssp = np.loadtxt(f"{fpath}/HalfSpace/RCVR{stn}_grav.h.dat")
        ax[irow,0].plot(qssp[:,0], qssp[:,index], color='blue', linewidth=lw, label='QSSP (Zhang)')

        ax[irow,0].set_xlim([-10, 50])
        ax[irow,0].set_ylabel('Gravity Strain (x 1e-12)')

    for irow in range(3):
        stn = irow + 4
        anal = np.loadtxt(f"{fpath}/Analytical/grstrain.RCVR{stn}.TOHOKU.dat")
        ax[irow,1].plot(anal[:,0], anal[:,index], color=analytical_color, linewidth=lw, label='Analytical Half Space')

        # Halfspace data QSSP:
        qssp = np.loadtxt(f"{fpath}/HalfSpace/RCVR{stn}_grav.h.dat")
        ax[irow,1].plot(qssp[:,0], qssp[:,index], color='blue', linewidth=lw, label='QSSP (Zhang)')


        ax[irow,1].set_xlim([-20, 150])



    ax[2,0].set_xlabel('Time [s]')
    ax[2,1].set_xlabel('Time [s]')

    return fig, ax










