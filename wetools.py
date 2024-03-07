import matplotlib.pyplot as plt
import numpy as np
import obspy
from obspy.core.trace import Trace
from obspy.core.stream import Stream


def moving_avg(f, time=[], half_window=1000, convert_t=False):
    # =====================================================================================================================================
    # DESCRIPTION:
    # Function produces a moving-average time series of user-inputted time-series with window-size 2*spacing+1
    # In the case that u = u(t), user may want accompanying t time-series to be calculated
    # INPUTS:
    #    f [1D array]             - number of seconds
    #    time (opt) [int]         - coarse-grain scale
    #    time (opt) [int, array]  - default set to 'No' indicates time array doesnt need to be calculated; otherwise 1D time array/time-series
    # OUTPUTS:
    #    cgu [1D array]           - coarse-grained time-series
    #    cg_time (opt) [1D array] - accompanying time array for time-series
    # =====================================================================================================================================


    # Initialise mov_avg array:
    avg = np.zeros(int(len(f) - 2 * half_window))

    # Moving-average calculations cuts off the first and last number of elements equal to the value of 'half_window'
    for ts_index in range(half_window, int(len(f) - half_window)):
        avg[ts_index - half_window]= 1 / (half_window * 2) * np.sum(f[ts_index - half_window:ts_index + half_window])

    # If convering a time series eg u(t), may wish to slice time array to same size:
    if convert_t==True:
        try:
            t_avg = time[half_window:int(len(f) - half_window)]
            return t_avg, avg
        except time == []:
            print('Error: No time array inputted')

    else:
        return avg



def slice_by_time(t, d, t_low, t_high):
    # Get indices of values to slice by:
    # Lower bound:
    lb = np.where(t>=t_low)
    lb = lb[0][0]
    # Upper bound:
    ub = np.where(t<=t_high)
    ub = ub[0][-1]
    return t[lb:ub+1], d[lb:ub+1]


def normalise(x):
    return x/np.amax(np.abs(x))

def normalise_after_time(t, x, min_time):
    m = np.amax(np.abs(x[t>min_time]))
    return x/m


# Converts obspy trace into x and y for mpl plotting
def obspy_gen_mpl(tr):
    x = np.linspace(0, tr.stats.npts*tr.stats.delta,  tr.stats.npts)
    y = tr.data
    return x,y


def plot_obspy_MPL(tr):
    x, y = obspy_gen_mpl(tr)

    fig, ax = plt.subplots()
    ax.plot(x,y)
    return fig


def calc_offset(src_lat, src_lon, stn_lat, stn_lon, geoco=1):
    # Using formula:
    # \Delta = acos(sin(src_lon)*sin(stn_lon) + cos(src_lon)*cos(stn_lon)*con(abs(src_lat - stn_lat)))
    pi = np.pi
    theta_src = (pi / 2) - np.arctan(geoco * np.tan(src_lat * pi / 180))  # Source colatitude
    theta_stn = (pi / 2) - np.arctan(geoco * np.tan(stn_lat * pi / 180))  # Station colatitude

    phi_src = src_lon * pi / 180  # Source longitude
    phi_stn = stn_lon * pi / 180  # Station longitude

    # Calculate epicentral distance $\Theta$ (scalar):
    a = np.cos(theta_stn) * np.cos(theta_src)
    b = np.sin(theta_stn) * np.sin(theta_src) * np.cos(phi_stn - phi_src)
    offset = np.arccos(a + b) * 180 / pi  # In degrees

    return offset



def save_figs_pdf(figs, pdf_name):
    import matplotlib.backends.backend_pdf
    pdf = matplotlib.backends.backend_pdf.PdfPages(pdf_name)
    for fig in figs:
        pdf.savefig(fig)
    pdf.close()




def rotate_stream(stream,  method, src, stn, geoco=1, overwrite_channel_names=False, invert_p=True, invert_t=True):
# Following the method used by JT in NMSYNG (DT98 eqn 10.3)
    if len(stream)==1:
        # Possible that only Z in stream!
        if 'Z' in stream[0].stats.channel:
            print(f"Only {stream[0].stats.channel} for {stream[0].id}")
            print(f"NO ROTATION PERFORMED")
            return stream
        else:
            raise ValueError("Only one trace in Stream and it isnt Z!")
    elif len(stream) > 3:
        print(f"Too many channels...")
        return 1


    # Get station coordinates:
    lat_stn = stn[0]
    lon_stn = stn[1]

    pi = np.pi
    lat_src = src[0]
    lon_src = src[1]


    # Convert geographic decimal --> geocentric radian values (if geoco==1 then geographic == geocentric):
    # Note here that theta is the CO-latitude
    theta_src = (pi/2) - np.arctan(geoco * np.tan( lat_src * pi/180 ))  # Source colatitude
    theta_stn = (pi/2) - np.arctan(geoco * np.tan( lat_stn * pi/180 ))  # Station colatitude

    phi_src   = lon_src*pi/180                                               # Source longitude
    phi_stn   = lon_stn*pi/180                                               # Station longitude

    # Calculate epicentral distance $\Theta$ (scalar):
    dist = np.arccos(np.cos(theta_stn)*np.cos(theta_src) + np.sin(theta_stn)*np.sin(theta_src)*np.cos(phi_stn - phi_src))

    rot1 = (1/np.sin(dist)) * \
           (np.sin(theta_stn)*np.cos(theta_src)  -  np.cos(theta_stn)*np.sin(theta_src)*np.cos(phi_stn - phi_src))

    rot2 = (1/np.sin(dist)) * (np.sin(theta_src)*np.sin(phi_stn - phi_src))

    # Conversion from RTP --> ZNE (where R=Z, 2D rotation) appears to use the following matrix:
    #   [N, E]' = [-rot1, -rot2; -rot2, rot1][T, P]' where T and P are theta, Phi
    #   Below we shall name the rotation matrix Q:
    # Hence to get the T and P matrix we should be multiplying [N,E] by the inverse of Q:
    Q    = np.array([[-rot1, -rot2], [rot2, -rot1]])
    Qinv = np.linalg.inv(Q)


    if method == "NE->RT":
        N = stream.select(component="N")[0].data
        E = stream.select(component="E")[0].data
        data_NE = np.array([N,E])
        data_TP = np.matmul(Qinv, data_NE)

        # Now writing back to stream:
        old_chls = ["N", "E"]
        new_chls = ["R", "T"]
        for i in range(2):
            if np.logical_and(new_chls[i] == "T", invert_p==True):
                data_TP[i, :] = data_TP[i,:]*(-1)
            if np.logical_and(new_chls[i] == "R", invert_t==True):
                data_TP[i, :] = data_TP[i,:]*(-1)

            stream.select(component=old_chls[i])[0].data = data_TP[i,:]
            if overwrite_channel_names:
                old_chl_name = stream.select(component=old_chls[i])[0].stats.channel

                if old_chl_name.find('E') != -1:
                    # Channel has E in it:
                    index = old_chl_name.find('E')
                    new_name = old_chl_name[:index] + 'T' + old_chl_name[index+1:]

                if old_chl_name.find('N') != -1:
                    # Channel has N in it:
                    index = old_chl_name.find('N')
                    new_name = old_chl_name[:index] + 'R' + old_chl_name[index+1:]

                stream.select(component=old_chls[i])[0].stats.channel = new_name

        return 0
    else:
        raise ValueError("Currently method must be NE->RT")





def rotate_trace(traceN, traceE, method, src, stn, geoco=1):
    # Following the method used by JT in NMSYNG (DT98 eqn 10.3)
    # Get station coordinates:
    lat_stn = stn[0]
    lon_stn = stn[1]

    pi = np.pi
    lat_src = src[0]
    lon_src = src[1]


    # Convert geographic decimal --> geocentric radian values (if geoco==1 then geographic == geocentric):
    # Note here that theta is the CO-latitude
    theta_src = (pi/2) - np.arctan(geoco * np.tan( lat_src * pi/180 ))  # Source colatitude
    theta_stn = (pi/2) - np.arctan(geoco * np.tan( lat_stn * pi/180 ))  # Station colatitude

    phi_src   = lon_src*pi/180                                               # Source longitude
    phi_stn   = lon_stn*pi/180                                               # Station longitude

    # Calculate epicentral distance $\Theta$ (scalar):
    dist = np.arccos(np.cos(theta_stn)*np.cos(theta_src) + np.sin(theta_stn)*np.sin(theta_src)*np.cos(phi_stn - phi_src))

    rot1 = (1/np.sin(dist)) * \
           (np.sin(theta_stn)*np.cos(theta_src)  -  np.cos(theta_stn)*np.sin(theta_src)*np.cos(phi_stn - phi_src))

    rot2 = (1/np.sin(dist)) * (np.sin(theta_src)*np.sin(phi_stn - phi_src))

    # Conversion from RTP --> ZNE (where R=Z, 2D rotation) appears to use the following matrix:
    #   [N, E]' = [-rot1, -rot2; -rot2, rot1][T, P]' where T and P are theta, Phi
    #   Below we shall name the rotation matrix Q:
    # Hence to get the T and P matrix we should be multiplying [N,E] by the inverse of Q:
    Q    = np.array([[-rot1, -rot2], [rot2, -rot1]])
    Qinv = np.linalg.inv(Q)


    if method == "NE->TP":
        data_NE = np.array([traceN, traceE])
        data_TP = np.matmul(Qinv, data_NE)

        # Now writing back to stream:
        old_chls = ["N", "E"]
        new_chls = ["T", "P"]
        for i in range(2):
            if new_chls[i] == "P":
                data_TP[i, :] = data_TP[i,:]*(-1)

    else:
        raise ValueError("Currently method must be NE->TP")

    return data_TP[0,:], data_TP[1,:]




def add_day_hour_lines(ax):
    # Adds vertical lines for each hour for a 0 to 86400 second x axis
    for i in range(4):
        ax.axvline(i*3600*6, color='r', linestyle='--', alpha=0.5)
        for j in range(6):
            ax.axvline((i * 3600*6) + (j*3600) , color='b', linestyle='--', alpha=0.2)



def bandpass_synth(d, fmin, fmax):
    # d should be a 2d array: storing time and timeseries arrays
    # Sort input
    d = np.array(d)
    shape = np.shape(d)
    t1 = shape[0]==2
    t2 = shape[1]==2
    assert(np.logical_or(t1, t2))
    if t1:
        d = d.T

    # Init. trace
    t  = Trace()
    dt   = d[1,0] - d[0,0]
    t.data        =  d[:,1]
    t.stats.delta =  dt

    # Bandpass
    t = t.filter(type='bandpass', freqmin=fmin, freqmax=fmax)

    return t




def STF_convolve(data_dt, data, stf_t, stf, timeshift=0.):

    # Ensure resampling of the STF to the data timestep:
    tr = obspy.Trace()
    tr.data = stf
    tr.stats.delta = np.mean(stf_t[1:] - stf_t[:-1])                      # Mean timestep - assumes constant timestepping.
    tr.resample(1/data_dt)  # resample at data freq

    # Resampled value:
    stf_rs = tr.data
    stf_rs_t = np.linspace(stf_t[0], stf_t[-1], tr.stats.npts)

    conv       = np.convolve(data, stf_rs, mode="same")                      # Convolve signal
    time_conv  = np.arange(len(conv))*data_dt  + timeshift                # Generate corresponding convolution time array

    output     = np.transpose(np.array([time_conv, conv]))                # Collate data into 2D array
    #output     = output[int(np.floor(len(stf_rs_t)/2)):, :]                  # Slice lower end of array
    #output     = output[:-len(stf_rs_t), :]                                  # Slice upper end of array

    return output





def gauss_STF_convolve(time, data, half_duration, alpha=1.628):
    dt        = time[1] - time[0]                                        # Timestep (assumes a regularly spaced dt)

    stf_t, stf = gen_gaussian(dt, half_duration, alpha)


    conv      = np.convolve(data, stf, mode="full")                      # Convolve signal
    time_conv = np.arange(len(conv))*dt - 1.5*half_duration              # Generate corresponding convolution time array

    output    = np.transpose(np.array([time_conv, conv]))                # Collate data into 2D array
    output    = output[int(np.floor(len(stf_t)/2)):, :]                  # Slice lower end of array
    output    = output[:-len(stf_t), :]                                  # Slice upper end of array

    return output



def bp_filter_ty(t, y, freqmin, freqmax, filter='bandpass', corners=4):
    # Filters a time series using Obspy
    dt = np.mean(t[1:] - t[:-1])

    tr = obspy.Trace()
    tr.data = y
    tr.stats.delta = dt
    T = tr.filter(type=filter, freqmin=freqmin, freqmax=freqmax, corners=corners, zerophase=True)

    t_out, y_out = obspy_gen_mpl(T)
    t_out += t[0]

    return t_out, y_out


def normalise_by_area(x,y):
    return y/np.trapz(y, x)


def gen_gaussian(dt, hdur, alpha=1.628):
    stf_t     = np.arange(-1.5*hdur, 1.5*hdur+dt, dt)  # Create time array for STF
    fact      = alpha/(((np.pi)**0.5)*hdur)                     # Gaussian pre-factor
    stf       = fact*np.exp(-(alpha*stf_t/hdur)**2)             # Gaussian STF
    return stf_t, stf


def pad_by_time(side, timeval, time, f, pad_time_bool=True):
    # Pad time series f (1D) with a certain number of seconds of zeros:

    dt = np.mean(time[1:] - time[:-1])

    pad_len = int(np.ceil(timeval/dt))
    print(f'Pad size: {pad_len}')
    pad_lower = np.zeros(pad_len)
    pad_upper = np.zeros(pad_len)

    tpad = (1+ np.arange(pad_len)) *dt
    tpad_lower = time[0] - tpad[::-1]
    tpad_upper = tpad + time[-1]

    if side == 'upper':
        pad_lower = []
        tpad_lower = []
    elif side == 'lower':
        pad_upper = []
        tpad_upper = []
    elif side == 'both':
        pass
    else:
        raise ValueError("side must be 'upper', 'lower' or 'both' ")

    # Pad the data
    f_new    = np.concatenate((pad_lower, f, pad_upper))
    if pad_time_bool:
        t = np.concatenate((tpad_lower, time, tpad_upper))
        return t, f_new, pad_len
    else:
        return f_new, pad_len




def data_to_trace(dt, d):
    tr = obspy.Trace()
    tr.data = d
    tr.stats.delta = dt
    return tr






def ms2_to_gal(x):
    return x/0.01

def ms2_to_milligal(x):
    return (x/0.01)*1e3

def ms2_to_microgal(x):
    return (x/0.01)*1e6

def ms2_to_nanogal(x):
    return (x/0.01)*1e9




def bp_filter(t, d, fmin, fmax, zerophase, corners):
    # Convert to stream:
    dt = np.mean(t[1:] - t[:-1])
    trace = data_to_trace(dt, d)

    trace = trace.filter(type='bandpass', freqmin=fmin, freqmax=fmax, zerophase=zerophase, corners=corners)

    # Convert back to x, y:
    x, y = obspy_gen_mpl(trace)
    # Update starttime to be consistent:
    x += t[0]

    return x,y





def filter_timeseries(t, d, type, freq, zerophase, corners):
    # Convert to stream:
    dt = np.mean(t[1:] - t[:-1])
    trace = data_to_trace(dt, d)

    if type == 'bandpass':
        trace = trace.filter(type='bandpass', freqmin=freq[0],
                             freqmax=freq[1], zerophase=zerophase,
                             corners=corners)
    elif type == 'highpass':
        trace = trace.filter(type='highpass', freq=freq,
                             zerophase=zerophase, corners=corners)
    elif type == 'lowpass':
        trace = trace.filter(type='lowpass', freq=freq,
                             zerophase=zerophase, corners=corners)
    else:
        raise ValueError('Type must be bandpass, highpass or lowpass')

    # Convert back to x, y:
    x, y = obspy_gen_mpl(trace)
    # Update starttime to be consistent:
    x += t[0]

    return x,y



# Converts obspy trace into x and y for mpl plotting
def obspy_gen_mpl_withtime(tr, starttime=0):
    x = np.linspace(0, tr.stats.npts*tr.stats.delta,  tr.stats.npts)
    y = tr.data

    if starttime != 0:
        x = x + starttime
    return x,y