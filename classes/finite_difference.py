from abc import ABC, abstractmethod
import os
import numpy as np
import math
from wetools.geodetics import delaz

class FD():
    def __init__(self, fdir='./', from_readme=False):
        self.fdir     = fdir
        self.stations = []
        self.nstns    = 0

        self.from_readme = from_readme

        self.readme_file   = self.fdir + '/README'
        self.spfm_stn_file = self.fdir + '/STATIONS'

        if self.from_readme:
            self.setup_from_readme()

    def add_station(self, stn):
        self.stations.append(stn)
        self.nstns += 1

    def set_stencil(self, stencil):
        self.stencil = stencil

    def write_specfem_stations(self):
        # Create directory if not present
        self._create_fdir()

        f = open(self.spfm_stn_file, "w")

        for istn in self.stations:
            # Loop stencil stations for this station
            for iss in range(istn.nst_stns):
                f.write(f"{istn.snms[iss]}  {istn.network}  {istn.slat[iss]:.6f}  {istn.slon[iss]:.6f}  0.0   0.0\n")
        f.close()

    def _create_fdir(self):
        # Creates file directory if not existing:
        if os.path.isdir(self.fdir):
            pass
        else:
            print(f'Creating directory: {self.fdir}')
            os.mkdir(self.fdir)


    def write_readme(self, **kwargs):

        defaultKwargs = {'fpath'   : self.readme_file}
        kwargs = defaultKwargs | kwargs

        f = open(kwargs['fpath'], "w")

        f.write(f"Number of stations: {self.nstns}\n")

        for istn in self.stations:
            istn._write_stn_to_readme(f)


    def setup_from_readme(self):
        f = open(self.readme_file, "r")
        l = f.readlines()

        nparams = 14

        # First line has number of stations:
        for istn in range(int(l[0].split()[-1])):

            ispar       = istn*nparams
            stnnme      =    str(l[ispar + 1].split()[-1])
            net         =    str(l[ispar + 2].split()[-1])
            clat        =  float(l[ispar + 3].split()[-1])
            clon        =  float(l[ispar + 4].split()[-1])
            dlat        =  float(l[ispar + 5].split()[-1])
            dlon        =  float(l[ispar + 6].split()[-1])
            rad         =  float(l[ispar + 7].split()[-1])
            geoco       =  float(l[ispar + 8].split()[-1])
            stencil     =   int(l[ispar + 10].split()[-2])
            order       =   int(l[ispar + 11].split()[-1])
            time_offset = float(l[ispar + 12].split()[-1])
            time_cutoff = float(l[ispar + 13].split()[-1])

            self.add_station(FdStation(stnname=stnnme,
                                       network=net,
                                       lat=clat,
                                       lon=clon,
                                       dlat=dlat,
                                       dlon=dlon,
                                       geoco=geoco,
                                       stencil=stencil,
                                       radius=rad,
                                       order=order,
                                       time_offset=time_offset,
                                       time_cutoff=time_cutoff))

    def load_specfem_data(self,types,dpath, acc_mult=1):
        for istn in self.stations:
            istn.load_specfem_data(types=types, dpath=dpath, acc_mult=acc_mult)




    def strain_from_acc(self, method, dx=None, code='specfem'):
        # Computes surface second gradient of gpot based on gacc (grad gpot):
        # dx should be in SI units (cartesian, not lat or lon)
        order = 1

        for istn in self.stations:

             if code =='specfem':
                d = istn.data
             elif code == 'qssp':
                d = istn.qssp
             else:
                 raise ValueError()

             for chl in ['NN', 'EE', 'ZN', 'EN', 'NE', 'ZE']:
                # First letter is data polarisation
                # Second letter is direction of dx
                c1, c2, nmp2, nmp1, nmm1, nmm2, ddx = \
                    self._get_fd_surface_station_strings(istn, chl, dx)

                datchl = f'{c1}C.PGRAV'
                outchl = f'{chl}.STRAIN'

                if method == 'centred_FD3':
                    d[outchl] = self.compute_derivative(order, method, dx=ddx,
                                                                 dp1=d[nmp1][datchl],
                                                                 dm1=d[nmm1][datchl])
                elif method == 'centred_FD5':
                    d[outchl] = self.compute_derivative(order, method, dx=ddx,
                                                                 dp1=d[nmp1][datchl],
                                                                 dp2=d[nmp2][datchl],
                                                                 dm1=d[nmm1][datchl],
                                                                 dm2=d[nmm2][datchl])
                else:
                    raise ValueError()

             if code =='specfem':
                 avg = (d[f'EN.STRAIN'] + d[f'NE.STRAIN'])/2
                 d[f'EN.STRAIN'] = avg
                 d[f'NE.STRAIN'] = avg





    def acc_from_gpot(self, method, dx=None):
        # Computes gradient of gpot based on gpot:
        # dx should be in SI units (cartesian, not lat or lon)
        order = 1

        for istn in self.stations:
            # We can compute NN and EE directly as 2nd derivatives
             for chl in ['N', 'E']:
                cc = istn.name + '_' + str(istn.cstn)
                datchl = f'G'
                outchl = f'{chl}C.PGRAV'

                c1, c2, nmp2, nmp1, nmm1, nmm2, ddx = \
                    self._get_fd_surface_station_strings(istn, chl='-'+chl, dx=dx)

                if method == 'centred_FD3':
                    istn.data[outchl] = self.compute_derivative(order, method, dx=ddx,
                                                                dc =istn.data[cc][datchl],
                                                                dp1=istn.data[nmp1][datchl],
                                                                dm1=istn.data[nmm1][datchl])
                elif method == 'centred_FD5':
                    istn.data[outchl] = self.compute_derivative(order, method, dx=ddx,
                                                                dc=istn.data[cc][datchl],
                                                                dp1=istn.data[nmp1][datchl],
                                                                dp2=istn.data[nmp2][datchl],
                                                                dm1=istn.data[nmm1][datchl],
                                                                dm2=istn.data[nmm2][datchl])
                else:
                    raise ValueError()







    def strain_from_gpot(self, method, dx=None):
        # Computes surface second gradient of gpot based on gpot:
        # dx should be in SI units (cartesian, not lat or lon)
        order = 2

        for istn in self.stations:
            # We can compute NN and EE directly as 2nd derivatives
             for chl in ['NN', 'EE']:
                # First letter is data polarisation
                # Second letter is direction of dx
                c1, c2, nmp2, nmp1, nmm1, nmm2, ddx = \
                    self._get_fd_surface_station_strings(istn, chl, dx)
                cc = istn.name + '_' + str(istn.cstn)

                datchl = f'G'
                outchl = f'{chl}.STRAIN'

                if method == 'centred_FD3':
                    istn.data[outchl] = self.compute_derivative(order, method, dx=ddx,
                                                                dc =istn.data[cc][datchl],
                                                                dp1=istn.data[nmp1][datchl],
                                                                dm1=istn.data[nmm1][datchl])
                elif method == 'centred_FD5':
                    istn.data[outchl] = self.compute_derivative(order, method, dx=ddx,
                                                                dc=istn.data[cc][datchl],
                                                                dp1=istn.data[nmp1][datchl],
                                                                dp2=istn.data[nmp2][datchl],
                                                                dm1=istn.data[nmm1][datchl],
                                                                dm2=istn.data[nmm2][datchl])
                else:
                    raise ValueError()



    def _get_fd_surface_station_strings(self, istn, chl, dx):
        # Determines the station codes for FD stations for computing a gradient
        # based on the channels (N,E,Z)
        c1 = chl[0];
        c2 = chl[1]

        if c2 == 'N':
            nmp1 = istn.name + '_' + str(int(istn.cstn + 10)).zfill(2)
            nmp2 = istn.name + '_' + str(int(istn.cstn + 20)).zfill(2)
            nmm1 = istn.name + '_' + str(int(istn.cstn - 10)).zfill(2)
            nmm2 = istn.name + '_' + str(int(istn.cstn - 20)).zfill(2)
            if dx == None:
                ddx = istn.dx_from_dlat()
            else:
                ddx = dx
        elif c2 == 'E':
            nmp1 = istn.name + '_' + str(int(istn.cstn + 1))
            nmp2 = istn.name + '_' + str(int(istn.cstn + 2))
            nmm1 = istn.name + '_' + str(int(istn.cstn - 1))
            nmm2 = istn.name + '_' + str(int(istn.cstn - 2))
            if dx == None:
                ddx = istn.dx_from_dlon()
            else:
                ddx = dx
        else:
            raise ValueError(f'C2 can only be North or East and is {c2}')

        return c1, c2, nmp2, nmp1, nmm1, nmm2, ddx


    def load_qssp_grav_acc(self,fpath):

        stnctr = 1
        for istn in self.stations:

            istn._check_qssp_exists()
            q = istn.qssp

            mult = 1  # They output gravity acceleration not grad phi hence the minus sign
            for chl in 'zne':
                # Columns are: time, centre, west, east, south, north
                d = np.loadtxt(fname=f'{fpath}/RCVR{stnctr}_grav_{chl}.dat', skiprows=1)

                timemask = d[:,0] < istn.time_cutoff
                q['time'] = d[timemask,0]

                chllbl = chl.upper() + 'C.PGRAV'
                q[istn.name + '_33'][chllbl] = mult*d[timemask, 1]  # Centre
                q[istn.name + '_32'][chllbl] = mult*d[timemask, 2]  # West
                q[istn.name + '_34'][chllbl] = mult*d[timemask, 3]  # East
                q[istn.name + '_23'][chllbl] = mult*d[timemask, 4]  # South
                q[istn.name + '_43'][chllbl] = mult*d[timemask, 5]  # North

            stnctr+=1

    def load_qssp_strain(self,fpath, coord_sys):
        stnctr = 1
        for istn in self.stations:

            istn._check_qssp_exists()
            q = istn.qssp

            if coord_sys=='RTZ':
                cchls = ['strain_time', '++', 'ZZ', 'XX', 'XX', 'RZ', 'TZ']
            elif coord_sys=='ENZ':
                cchls = ['strain_time', 'EE', 'NN', 'ZZ', 'EN', 'NE', 'ZE', 'ZN']
            else:
                raise ValueError(f'coord_sys can only be RTZ or ENZ but is {coord_sys}')

            d = np.loadtxt(fname=f'{fpath}/RCVR{stnctr}_grav.h.{coord_sys}.dat', skiprows=1)
            q['strain_time'] = d[:, 0]  # time
            for i in range(1,len(cchls)):
                q[f"{cchls[i]}.STRAIN"]    =  d[:, i]

            stnctr += 1



    def compute_derivative(self, order, method, dx, **kwargs):
        defaultKwargs = {'dc': None, 'dp2': None, 'dm2': None}
        kwargs = defaultKwargs | kwargs

        dc  = kwargs['dc']      # data at central point
        dp1 = kwargs['dp1']     # data at + delta x
        dm1 = kwargs['dm1']     # data at - delta x
        dp2 = kwargs['dp2']     # data at +  2 delta x
        dm2 = kwargs['dm2']     # data at -  2 delt

        # Centred FD:
        if method == 'centred_FD3':
            return self.compute_centred_FD3(order=order, p=dp1, m=dm1, dx=dx, c=dc)
        elif method == 'centred_FD5':
            return self.compute_centred_FD5(order=order, dx=dx, c=dc,
                                            p1=dp1, m1=dm1,
                                            p2=dp2, m2=dm2)
        else:
            raise ValueError(f"method {method} not available ")


    def compute_centred_FD3(self, order, p, m, dx, c=None):
        if order   == 1:
            # First order derivative
            deriv = (1/(2*dx)) * (p - m)
        elif order == 2:
            # Second order derivative
            if (c == None).any():
                raise ValueError('Data at central point is required.')
            deriv = (p + m - 2*c)/(dx*dx)
        return deriv


    def compute_centred_FD5(self, order, p1, p2, m1, m2, dx, c=None):
        # https://en.wikipedia.org/wiki/Five-point_stencil
        if order   == 1:
            # First order derivative
            deriv = (-p2 + 8*p1 - 8*m1 + m2)/(12*dx)
        elif order == 2:
            # Second order derivative
            if (c == None).any():
                raise ValueError('Data at central point is required.')
            deriv = (-p2 + 16*p1 - 30*c + 16*m1 - m2)/(12*dx*dx)
        elif order == 3:
            # Third order derivative
            deriv = (p2 - 2*p1 + 2*m1 - m2)/(2*dx*dx*dx)
        elif order == 4:
            # Fourth order derivative
            if (c == None).any():
                raise ValueError('Data at central point is required.')
            deriv = (p2 - 4*p1 + 6*c -4*m1 + m2)/(dx**4)
        else:
            raise ValueError(f'Order {order} isnt available. Only 1,2,3,4.')
        return deriv



    def rotate_strain(self, from_to, event_coords):

        if from_to == 'nez_to_rtz':

            zhang_rotate_strain(h_EE,
                                h_NN,
                                h_EN,
                                h_NE,
                                h_EZ,
                                h_NZ,
                                baz)

        else:
            raise ValueError('rotate_strain only nez_to_rtz method implemented')


class FdStation():

    def __init__(self, stnname, lat, lon, dlat, dlon, **kwargs):

        # Default kwargs:
        defaultKwargs = {'geoco'   : 0.9966471893352525,
                         'radius'  : 6371000.0,
                         'stencil' : 0,
                         'order'   : 0,
                         'network' : 'YY',
                         'time_offset': 0,
                         'time_cutoff': 1e30}

        kwargs = defaultKwargs | kwargs

        # Possible kwargs:
        self.geoco   = kwargs['geoco']
        self.radius  = kwargs['radius']
        self.stencil = kwargs['stencil']
        self.order   = kwargs['order']
        self.network = kwargs['network']
        self.time_offset = kwargs['time_offset']
        self.time_cutoff = kwargs['time_cutoff']

        self.name = stnname

        self.clat = lat
        self.clon = lon

        self.dlat = dlat
        self.dlon = dlon

        self.data = {}

        # Determine perfect sphere:
        if self.geoco == 1:
            self.perfect_sphere = True
        else:
            self.perfect_sphere = False

        if self.stencil > 0 and self.order > 0:
            # Have parsed stencil/order so create:
            if    self.stencil == 3:
                self.stenciltype  = Stencil3
                self.cstn = 11
            elif  self.stencil == 5:
                self.stenciltype  = Stencil5
                self.cstn = 33

            self.cstnname = self.name + '_' + str(self.cstn)

            self.stencil = self.stenciltype(order=self.order)

            # Apply the stencil:
            self.apply_stencil_to_cstn()

    def dx_from_dlat(self):
        # Converts latitude distance to metres
        if self.perfect_sphere:
            return self.spherical_to_cart_dis(angle=self.dlat, geoco=self.geoco)
        else:
            raise ValueError('Not implemented for non-spherical body')

    def dx_from_dlon(self):
        # Converts longitude distance to metres
        if self.perfect_sphere:
            return self.spherical_to_cart_dis(angle=self.dlon, geoco=self.geoco)
        else:
            raise ValueError('Not implemented for non-spherical body')


    def compute_epicentral_distance(self, srclat, srclon):
        # Compute the azimuthal surface distance for each station:
        self.epi_dist, azz, azz1 = delaz(eplat=srclat, eplong=srclon, stlat=self.clat, stlong=self.clon, geoco=self.geoco)
        self.epi_dist_km = self.spherical_to_cart_dis(self.epi_dist, geoco=self.geoco) / 1000 # km


    def spherical_to_cart_dis(self, angle, geoco=1):
        return angle*(np.pi*self.radius)/180


    def apply_stencil_to_cstn(self):

        self.nst_stns = 2*self.stencil.points - 1

        self.slat = []
        self.slon = []
        self.snms = []

        frommid = int((self.stencil.points -1)/2)
        for ii in range(-frommid, frommid+1):
            self.slat.append(np.around(self.clat + ii*self.dlat, 6))
            self.slon.append(np.around(self.clon, 6))

        for jj in range(-frommid, frommid+1):
            if jj != 0:
                self.slat.append(np.around(self.clat, 6))
                self.slon.append(np.around(self.clon + jj * self.dlon, 6))

        for kk in range(self.nst_stns):
            self.snms.append(f"{self.name}_{self.stencil.stencil[kk]}")


    def _write_stn_to_readme(self, f):

        f.write(f"Station                :        {self.name}\n")
        f.write(f"  - network            :        {self.network}\n")
        f.write(f"  - central latitude   :        {self.clat}\n")
        f.write(f"  - central longitude  :        {self.clon}\n")
        f.write(f"  - delta latitude     :        {self.dlat}\n")
        f.write(f"  - delta longitude    :        {self.dlon}\n")

        f.write(f"  - radius             :        {self.radius}\n")
        f.write(f"  - geoco              :        {self.geoco}\n")
        f.write(f"  - perfect sphere     :        {self.perfect_sphere}\n")
        f.write(f"  - stencil            :        {self.stencil.points} point\n")
        f.write(f"  - order              :        {self.order}\n")
        f.write(f"  - time offset        :        {self.time_offset}\n")
        f.write(f"  - time cutoff        :        {self.time_cutoff}\n\n")


    def load_specfem_data(self, types, dpath, acc_mult=1):

        for stnname in self.snms:
            stndata = {}

            loadpref = f"{dpath}/{stnname}.{self.network}"

            for type in types:
                loadstr = []
                labels = []

                if type == 'pgrav':
                    for chl in 'ZNE':
                        label = f"{chl}C.PGRAV"
                        loadstr.append(f'{loadpref}.BX{label}.sem.ascii')
                        labels.append(label)
                elif type == 'strain':
                    for chl in ['ZZ', 'NN', 'EE', 'NZ', 'NE', 'EZ']:
                        label = f"{chl}.STRAIN"
                        loadstr.append(f'{loadpref}.BX{label}.sem.ascii')
                        labels.append(label)
                elif type == 'gpot':
                        label = "G"
                        loadstr.append(f'{loadpref}.BX{label}.sem.ascii')
                        labels.append(label)
                else:
                    raise ValueError(f'{type} data not implemented yet')

                # Load data:
                for ic in range(len(loadstr)):
                    d = np.loadtxt(loadstr[ic])
                    d[:, 0] += self.time_offset

                    timemask = d[:,0]  < self.time_cutoff

                    if type == 'pgrav':
                        d[:,1]*=acc_mult

                    stndata[labels[ic]] = d[timemask,1]


            self.data['time']  = d[timemask,0]
            self.data[stnname] = stndata

    def _check_qssp_exists(self):
        # See if dict exists
        try:
            self.qssp
        except:
            self.qssp = {}
            for sstn in ['_33','_32','_34','_23','_43']:
                self.qssp[self.name + sstn] = {}



class Stencil(ABC):
    @abstractmethod
    def __init__(self, order):
        self.order  = order

    def _setup_stencil(self):
        self.stencil      = []

        # Set up stencil:
        for ilat in range(self.points):
                self.stencil.append(f'{self.start + ilat}{self.midpoint}')
        for ilon in range(self.points):
            if ilon + self.start != self.midpoint:
                self.stencil.append(f'{self.midpoint}{self.start + ilon}')


class Stencil3(Stencil):

    def __init__(self, order):
        super().__init__(order)

        self.points   = 3
        self.midpoint = 1
        self.start    = 0
        self._setup_stencil()


class Stencil5(Stencil):

    def __init__(self, order):
        super().__init__(order)

        self.points   = 5
        self.midpoint = 3
        self.start    = 1
        self._setup_stencil()






def zhang_rotate_strain(h_EE, h_NN, h_EN, h_NE, h_EZ, h_NZ, baz):
    """
    rotate hxyz to hrtz according to back-azimuth by Shenjian Zhang
    :param hxyz: hxyz: hxx, hyy, hxy, hyx, hxz, hyz
    :param baz: backazimuth
    :return: hrtz: hrr, htt, hrt, htr, hrz, htz
    """
    # gravity strain: Nx6
    hxyz = np.stack((h_EE, h_NN, h_EN, h_NE, h_EZ, h_NZ), axis=1)

    alpha = 270-baz    # in degree
    ca = math.cos(alpha/180*math.pi)
    sa = math.sin(alpha/180*math.pi)
    hrr = ca*ca*hxyz[:,0] + sa*sa*hxyz[:,1] + 2*ca*sa*hxyz[:,2]
    htt = sa*sa*hxyz[:,0] + ca*ca*hxyz[:,1] - 2*ca*sa*hxyz[:,2]
    hrt = -ca*sa*hxyz[:,0] + ca*sa*hxyz[:,1] + (ca*ca-sa*sa)*hxyz[:,2]
    htr = hrt
    hrz = ca*hxyz[:,4] + sa*hxyz[:,5]
    htz = -sa*hxyz[:,4] + ca*hxyz[:,5]

    hrtz = np.stack((hrr,htt,hrt,htr,hrz,htz), axis=1)
    return hrtz