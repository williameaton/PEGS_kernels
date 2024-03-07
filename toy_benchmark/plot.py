import numpy as np 
import matplotlib.pyplot as plt 

net = 'Y5'
stn = 'Z1'
chls = ['Z', 'N', 'E', 'G']

nochls= len(chls)
fig, ax = plt.subplots(nochls+1,2, figsize=(12, 7.5), sharex=True)
fig.set_tight_layout(True)
chlctr = 0 

plot_stored_forward        = True
plot_reconstructed_forward = True

plot_stable_adjoint     = True
plot_calculated_adjoint = True 
# Reconstructed forward traces: 
for chl in chls: 

    # Load reconstructed simulation data:
    if plot_reconstructed_forward: 
        if chl=='G': 
            d = np.loadtxt(f"./STORE_RECONSTRUCTED_G/{stn}.{net}.MX{chl}.sem.ascii")
            # Stored copy for estimating error
            gt = d[:, 0]
            gs = d[:,1]
            lenG = len(gs)

        else: 
            d = np.loadtxt(f"./STORE_RECONSTRUCTED_G/{stn}.{net}.MX{chl}D.sem.ascii")
        ax[chlctr,0].plot(d[:,0], d[:,1], 'k-', linewidth=0.4)


    if plot_stored_forward: 
        # load the stored forward traces from the original run: 
        d = np.loadtxt(f"./SEM/{stn}.{net}.MX{chl}.adj")
        if chl=='G': 
            ax[chlctr,0].plot(d[::-1,0], d[:,1], 'r:', linewidth=0.6)
            stt = d[:, 0]
            sst = d[::-1, 1]
        else: 
            ax[chlctr,0].plot(d[::-1,0], d[:,1], 'r:', linewidth=0.6)


    ax[chlctr,0].set_title(f'Channel: {chl}')
    ax[0,0].legend(['Reconstructed forward', 'Original forward'])

    if chl =='G': 
        misfit = np.abs((gs[:lenG]- sst[:lenG])/sst[:lenG])
    

        ax[-1,0].plot(gt[:lenG], misfit , 'kx')

        ax[-1,0].set_title('log abs fractional (max is 1 = 100 %) difference of G traces')
        ax[-1,0].set_yscale('log')
        ax[-1,0].set_ylim([1e-7, 1e-1])
        ax[-1,0].axhline(1e-2, color='k', alpha=0.3, linestyle='--', linewidth=0.5)
        ax[-1,0].axhline(1e-3, color='k', alpha=0.3, linestyle='--', linewidth=0.5)
        ax[-1,0].axhline(1e-4, color='k', alpha=0.3, linestyle='--', linewidth=0.5)
        ax[-1,0].axhline(1e-5, color='k', alpha=0.3, linestyle='--', linewidth=0.5)
        ax[-1,0].axhline(1e-6, color='k', alpha=0.3, linestyle='--', linewidth=0.5)

    # Load adjoint seismogram from current output run:
    lgnd = []
    lgndstr = []
    if plot_calculated_adjoint:
        if chl=='G': 
            d = np.loadtxt(f"./OUTPUT_FILES/{stn}.{net}.MX{chl}.sem.ascii")
            # Stored copy for estimating error
            gt = d[:, 0]
            gs = d[:,1]
        else: 
            d = np.loadtxt(f"./OUTPUT_FILES/{stn}.{net}.MX{chl}D.sem.ascii")
        lc, = ax[chlctr,1].plot(d[:,0], d[:,1], 'k', linewidth=0.6)
        lgnd.append(lc)
        lgndstr.append('Current adjoint')

    if plot_stable_adjoint: 
        # load the stored forward traces from the original run: 
        if chl=='G': 
            d = np.loadtxt(f"./STABLE_ADJOINT_TRACES/{stn}.{net}.MX{chl}.sem.ascii")
            stt = d[:, 0]
            sst = d[::-1, 1]
            lenG = len(stt)

        else: 
            d = np.loadtxt(f"./STABLE_ADJOINT_TRACES/{stn}.{net}.MX{chl}D.sem.ascii")
        lc, = ax[chlctr,1].plot(d[:,0], d[:,1], 'r--', linewidth=0.6)
        lgnd.append(lc)
        lgndstr.append('Stable adjoint')


    if chl =='G': 
        #misfit = np.abs((gs[:lenG]- sst[:lenG])/sst[:lenG])
        #ax[-1,1].plot(gt[:lenG], misfit , 'kx')

        ax[-1,1].set_title('log abs fractional (max is 1 = 100 %) difference of G traces')
        ax[-1,1].set_yscale('log')
        ax[-1,1].set_ylim([1e-7, 1e-1])
        ax[-1,1].axhline(1e-2, color='k', alpha=0.3, linestyle='--', linewidth=0.5)
        ax[-1,1].axhline(1e-3, color='k', alpha=0.3, linestyle='--', linewidth=0.5)
        ax[-1,1].axhline(1e-4, color='k', alpha=0.3, linestyle='--', linewidth=0.5)
        ax[-1,1].axhline(1e-5, color='k', alpha=0.3, linestyle='--', linewidth=0.5)
        ax[-1,1].axhline(1e-6, color='k', alpha=0.3, linestyle='--', linewidth=0.5)

    ax[chlctr,1].set_title(f'Channel: {chl}')
    ax[0,1].legend(lgnd, lgndstr)





    chlctr += 1 






plt.savefig(f"{stn}.{net}.pdf", format="pdf", bbox_inches="tight")
