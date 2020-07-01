import numpy as np


def isfloat(value):
    try:
        float(value)
        return True
    except ValueError:
        return False


class PincFloitIO:
    """
properties:
 - nx, ny, nz : number of grid cells in x,y,z-direction
 - dt_write : length of output interval in s
 - t_max : total time of simulation in s
 - xlim, ylim, zlim = array containing upper and lower domain bounds in x,y,z-direction in m
 - nvar : number of variables
 - nt : number of time steps
 - tt : array of points in time
 - xx, yy, zz: spatial coordinates in array form
    """

    def __init__(self, path, configpath, logpath=None, dellog=False, nvar=None, nframe=None):
        self.raw_data = np.fromfile(path, dtype='float32')

        with open(configpath) as file:
            self.config = file.readlines()

        for nn in range(0, len(self.config)):
            line = self.config[nn]
            if (line[0] != '!') and ('sizeX' in line):
                self.nx = \
                    [int(substring.replace(',', '')) for substring in line.split() if
                     isfloat(substring.replace(',', ''))][0]
            if (line[0] != '!') and ('sizeY' in line):
                self.ny = \
                    [int(substring.replace(',', '')) for substring in line.split() if
                     isfloat(substring.replace(',', ''))][0]
            if (line[0] != '!') and ('sizeZ' in line):
                self.nz = \
                    [int(substring.replace(',', '')) for substring in line.split() if
                     isfloat(substring.replace(',', ''))][0]
            if (line[0] != '!') and ('outputTimeDiff' in line):
                self.dt_write = [float(substring.replace(',', ''))
                                 for substring in line.split() if isfloat(substring.replace(',', ''))][0]
            if (line[0] != '!') and ('maxTime' in line):
                self.t_max = [float(substring.replace(',', '')) for substring in line.split() if
                              isfloat(substring.replace(',', ''))][0]
            if (line[0] != '!') and ('lx_dim' in line):
                self.xlim = [float(substring.replace(',', '')) for substring in line.split() if
                             isfloat(substring.replace(',', ''))]
            if (line[0] != '!') and ('ly_dim' in line):
                self.ylim = [float(substring.replace(',', '')) for substring in line.split() if
                             isfloat(substring.replace(',', ''))]
            if (line[0] != '!') and ('lz_dim' in line):
                self.zlim = [float(substring.replace(',', '')) for substring in line.split() if
                             isfloat(substring.replace(',', ''))]

        self.nvar = None
        if nvar and nframe:
            self.nvar = nvar
            self.t_max = self.dt_write * nframe
        else:
            for nn in range(1, 15):
                if self.raw_data.shape[0] / self.nx / self.ny / self.nz / nn == int(self.t_max / self.dt_write) + 1:
                    self.nvar = nn

        if not self.nvar:
            raise Exception('Could not determine number of variables. Please Check the input.')

        self.nt = int(np.round(self.t_max / self.dt_write + 1)) # only works  if run was not interrupted!
        try:
            self.data = np.reshape(self.raw_data, (self.nt, self.nvar, self.nz, self.ny, self.nx))
        except IOError:
            raise Exception('Shape of data is wrong. Please check input files and specify nvar and nframe correctly.')

        if not logpath:
            self.tt = np.arange(0, self.nt) * self.dt_write
        else:
            with open(logpath) as file:
                self.log = file.readlines()

            aux = []
            for nn in range(0, len(self.log)):
                line = self.log[nn]
                if (line[0] != '!') and ('at physical time' in line):
                    aux.append([float(substring.replace(',', '')) for substring in line.split() if
                                isfloat(substring.replace(',', ''))][0])
            self.tt = np.array(aux)

            if dellog:
                del self.log

        self.xx = np.linspace(self.xlim[0], self.xlim[1], self.nx)
        self.yy = np.linspace(self.ylim[0], self.ylim[1], self.ny)
        self.zz = np.linspace(self.zlim[0], self.zlim[1], self.nz)
        self.namelist = None

    def set_varnames(self, namelist):
        """
        Set variable names for variables in file from array-like pbject namelist

        :param namelist - list of variables / attributes to the object
        """
        if len(namelist) != self.nvar:
            raise Exception('Number of supplied variable names does not match number of variables!')

        setattr(self, 'namelist', namelist)
        for nn in range(0, len(namelist)):
            setattr(self, namelist[nn], self.data[:, nn])

    def ret(self, varname):  # return attribute from string of attribute name
        """
        returns variable with given name; this is intended to be used for calling an attribute with a string

        :param varname - name of variable [str]
        """
        return self.__getattribute__(varname)

    def get_energy(self, nt, theta_0=300, zz_0=40e3, NN_0=1e-2, gg=9.81, mlimits=None):
        """
        Get the IGW energy by projection onto IGW mode in constant stratification.

        :param nt:          frame to be analyzed
        :param theta_0:     reference potential temperature
        :param zz_0:        reference height
        :param NN_0:        reference stratification
        :param gg:          gravity
        :param mlimits:     array with limits in vertical wave number to split
                            each element needs two components: lower limit / upper limit: shape(n, 2)
        :return:            IGW energy field; array-like according to wave number limits
        """
        print('{0:2d} / {1:2d}'.format(nt + 1, self.nt))

        theta_bg = theta_0 * np.exp(NN_0 ** 2 / gg * (self.zz - zz_0))

        # get wave numbers
        mm_1d = np.fft.fftfreq(self.nz, d=np.diff(self.zz[:2])[0])
        kx_1d = np.fft.fftfreq(self.nx, d=np.diff(self.xx[:2])[0])

        # set limits in vertical wave number for integration
        if mlimits is None:
            mlimits = np.array([[np.min(mm_1d), np.max(mm_1d)]])

        kx, mm = np.meshgrid(kx_1d, mm_1d)
        omega_p = NN_0 * (kx ** 2 / (kx ** 2 + mm ** 2)) ** .5
        omega_m = -omega_p
        # omega_m[0, 0] = omega_m[0, 0] = np.nan

        # first transform the pressure to construct the modes
        pp = np.ones((self.nz, self.nx))

        norm = (kx ** 2 / omega_p ** 2 + mm ** 2 * (omega_p ** 2 + NN_0 ** 2) / (
                NN_0 ** 2 - omega_p ** 2) ** 2) ** .5 * np.abs(pp) * .5 ** .5

        # calculate modes
        modes_p = ((np.ones(pp.shape + (4,)).transpose(2, 0, 1) *
                    np.array([kx / omega_p,
                              np.zeros(kx.shape),
                              - mm * omega_p / (NN_0 ** 2 - omega_p ** 2),
                              1j * mm * NN_0 / (NN_0 ** 2 - omega_p ** 2)
                              ])) * pp / norm)

        modes_m = ((np.ones(pp.shape + (4,)).transpose(2, 0, 1) *
                    np.array([kx / omega_m,
                              np.zeros(kx.shape),
                              - mm * omega_m / (NN_0 ** 2 - omega_m ** 2),
                              1j * mm * NN_0 / (NN_0 ** 2 - omega_m ** 2)
                              ])) * pp / norm)

        # set coefficient to zero where the pressure is zero / modes are NaN
        for nn in range(0, 4):
            modes_p[np.where(np.isnan(modes_p))] = 0.0
            modes_m[np.where(np.isnan(modes_m))] = 0.0

        coeff_p = np.zeros(pp.shape).astype(complex)
        coeff_m = np.zeros(pp.shape).astype(complex)

        vars = ['uu', 'vv', 'ww', 'theta']
        for nn in range(0, len(vars)):
            var = vars[nn]
            if var == 'theta':
                # construct buoyancy from potential temperature
                transform = np.fft.fft2(
                    gg * (self.ret(var)[nt, :, 0, :].transpose(1, 0) / theta_bg).transpose(1, 0)) / NN_0
            else:
                transform = np.fft.fft2(self.ret(var)[nt, :, 0, :])

            coeff_m += .5 * transform * modes_m[nn].conj()
            coeff_p += .5 * transform * modes_p[nn].conj()

        # construct spectral wave field
        zz_hat = coeff_p * modes_p + coeff_m * modes_m

        # calculate energy content in IGW field
        energy = np.empty((len(mlimits), self.nz, self.nx))
        for nn in range(0, len(mlimits)):

            index = np.where((np.abs(mm_1d) <= mlimits[nn, 0]) | (np.abs(mm_1d) > mlimits[nn, 1]))
            temp = zz_hat.copy()

            temp[:, index] = 0.0
            zz = np.fft.ifft2(temp, axes=(1, 2))
            energy[nn] = .5 * (zz * zz.conj()).sum(axis=0).real

        return energy
