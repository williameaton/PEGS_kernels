import numpy as np
import matplotlib.pyplot as plt

fig, ax = plt.subplots(3)

ichl = 0
for chl in 'ZNE':
    stn = 'ULN'

    # Load the adjoint source:
    adj = np.loadtxt(f'{stn}_DATA/adjoint_source_{stn}_{chl}')
    adj_time = adj[0,:] - 105
    adj_src  = adj[1,:]

    # Now load original time series to get exact times needed:
    spfm = np.loadtxt(f'{stn}_DATA/plot_source_time_function.txt')


    ax[ichl].plot(adj_time, adj_src)

    #ax[ichl].plot(spfm[:,0], spfm[:,0]*0, 'x')


    newt   = spfm[:,0]
    newsrc = np.zeros(len(spfm))


    for i in range(1, len(adj_src)):
        tnew = newt[i]

        # Find the closest two indexes to this time, on the old time array
        mintime = adj_time[adj_time < tnew][-1]
        mindex  = int(np.where(adj_time==mintime)[0][0])

        maxtime = adj_time[adj_time > tnew][0]
        pindex  = int(np.where(adj_time==maxtime)[0][0])


        newval = adj_src[mindex] +  (tnew - mintime)*(adj_src[pindex]-adj_src[mindex])/(maxtime-mintime)

        newsrc[i] = newval



    # Source is negative derivative so flip sign:
    newsrc *=-1

    ax[ichl].plot(newt, newsrc)

    ax[ichl].axvline(newt[184])
    ax[ichl].axvline(newt[-185])



    np.savetxt(f'{stn}_DATA/{stn}_0.IC.MX{chl}.adj', X=np.array([newt, newsrc]).T)

    ichl += 1

plt.show()