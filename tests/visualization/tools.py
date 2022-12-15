"""
This Python module provides tools for importing and transforming model output
data from PincFlow.
"""

import numpy
import xarray


class ModelOutput():
    """Import and transform model output data from PincFlow."""

    def import_data(self, path):
        """Import and format data."""

        # Read input file.
        with open("".join((path, "input.f90"))) as file:
            lines = file.readlines()

        # Remove white spaces, empty lines, comments and namelist titles.
        lines = [entry.split("!", 1)[0].strip() for entry in lines
                if len(entry.strip()) > 1 and entry.strip()[0] not
                in ("!", "&")]

        # Get keys and values.
        keys = [entry.split("=", 1)[0].strip() for entry in lines]
        values = [entry.split("=", 1)[1].strip() for entry in lines]

        # Remove trailing commas.
        values = [entry.rstrip(",") if "," in entry else entry for entry
                in values]

        # Adjust boolean entries.
        values = [entry.replace(".true.", "True") if ".true." in entry
                else entry for entry in values]
        values = [entry.replace(".false.", "False") if ".false." in entry
                else entry for entry in values]

        # Evaluate entries.
        values = [eval(entry) if "=" not in entry else entry for entry
                in values]

        # Save to dictionary.
        input = {keys[nn]: values[nn] for nn in range(len(keys))}

        # Set attributes.
        self.nx = input["sizeX"]
        self.ny = input["sizeY"]
        self.nz = input["sizeZ"]
        self.lx = input["lx_dim"][1] - input["lx_dim"][0]
        self.ly = input["ly_dim"][1] - input["ly_dim"][0]
        self.lz = input["lz_dim"][1] - input["lz_dim"][0]
        self.n0 = input["N_BruntVaisala_dim"]
        self.f0 = input["f_Coriolis_dim"]
        self.u0 = input["backgroundFlow_dim"][0]
        self.v0 = input["backgroundFlow_dim"][1]
        self.w0 = input["backgroundFlow_dim"][2]
        self.topography = input["topography"]
        self.output_type = input["outputType"]
        self.dt = input["outputTimeDiff"]
        self.lt = input["maxTime"]
        self.npsi = numpy.sum(input["varOut"][:])

        # Set topography attributes.
        if self.topography:
            self.h0 = input["mountainHeight_dim"]
            self.l0 = input["mountainWidth_dim"]

        # Define grid spacing.
        self.dx = self.lx / self.nx
        self.dy = self.ly / self.ny
        self.dz = self.lz / self.nz

        # Define grid.
        self.xx = numpy.linspace(0.5 * self.dx, self.lx - 0.5 * self.dx,
                self.nx)
        self.yy = numpy.linspace(0.5 * self.dy, self.ly - 0.5 * self.dy,
                self.ny)
        self.zz = numpy.linspace(0.5 * self.dz, self.lz - 0.5 * self.dz,
                self.nz)
        self.zz, self.yy, self.xx = numpy.meshgrid(self.zz, self.yy, self.xx,
                indexing = "ij")

        # Import data.
        psi = numpy.fromfile("".join((path, "pf_all.dat")), dtype = "float32")
        self.nt = int(len(psi) / self.nx / self.ny / self.nz / self.npsi)
        self.psi = numpy.reshape(psi, (self.nt, self.npsi, self.nz, self.ny,
                self.nx))

        # Interpolate velocities to cell centers.
        if self.nx > 1:
            self.psi[:, 1, :, :, 1:] = 0.5 * (self.psi[:, 1, :, :, 1:]
                    + self.psi[:, 1, :, :, :(- 1)])
        if self.ny > 1:
            self.psi[:, 2, :, 1:, :] = 0.5 * (self.psi[:, 2, :, 1:, :]
                    + self.psi[:, 2, :, :(- 1), :])
        if self.nz > 1:
            self.psi[:, 3, 1:, :, :] = 0.5 * (self.psi[:, 3, 1:, :, :]
                    + self.psi[:, 3, :(- 1), :, :])

        # Compute output times.
        if self.output_type == "time":
            self.tt = numpy.zeros(self.nt)
            start = int(self.nt - self.lt / self.dt)
            if start >= 0:
                self.tt[start:] = numpy.arange(start, self.nt) * self.dt

        # Define topography.
        if self.topography:
            hh = numpy.fromfile("".join((path, "topography.dat")),
                    dtype = "float32")
            self.hh = numpy.reshape(hh, (self.ny, self.nx))

    def transform_data(self, interpolation = False):
        """Transform and interpolate data."""

        if self.topography:

            # Define Jacobian.
            jj = (self.lz - self.hh) / self.lz

            # Define metric tensor elements.
            gg = numpy.zeros((2, self.nz, self.ny, self.nx))
            gg[0, :, :, 1:(- 1)] = 0.5 * (self.hh[:, 2:]
                    - self.hh[:, :(- 2)]) / self.dx * (self.zz[:, :, 1:(- 1)]
                    - self.lz) / (self.lz - self.hh[:, 1:(- 1)])
            gg[1, :, 1:(- 1), :] = 0.5 * (self.hh[2:, :]
                    - self.hh[:(- 2), :]) / self.dy * (self.zz[:, 1:(- 1), :]
                    - self.lz) / (self.lz - self.hh[1:(- 1), :])

            # Define metric tensor elements at the horizontal boundaries.
            if self.nx > 1:
                gg[0, :, :, 0] = (0.5 * (self.hh[:, 1] - self.hh[:, - 1])
                        / self.dx * (self.zz[:, :, 0] - self.lz) / (self.lz
                        - self.hh[:, 0]))
                gg[0, :, :, - 1] = (0.5 * (self.hh[:, 0] - self.hh[:, - 2])
                        / self.dx * (self.zz[:, :, - 1] - self.lz) / (self.lz
                        - self.hh[:, - 1]))
            if self.ny > 1:
                gg[1, :, 0, :] = (0.5 * (self.hh[1, :] - self.hh[- 1, :])
                        / self.dy * (self.zz[:, 0, :] - self.lz) / (self.lz
                        - self.hh[0, :]))
                gg[1, :, - 1, :] = (0.5 * (self.hh[0, :] - self.hh[- 2, :])
                        / self.dy * (self.zz[:, - 1, :] - self.lz) / (self.lz
                        - self.hh[- 1, :]))

            # Compute Cartesian height.
            self.zc = jj * self.zz + self.hh

            # Compute Cartesian vertical wind.
            self.psi[:, 3, :, :, :] = (jj * self.psi[:, 3, :, :, :] - jj
                    * gg[0, :, :, :] * self.psi[:, 1, :, :, :] - jj
                    * gg[1, :, :, :] * self.psi[:, 2, :, :, :])

            # Interpolate to Cartesian grid (NumPy is way faster than SciPy).
            # This is not recommended for large arrays!
            if interpolation:
                psi = [[[[numpy.interp(self.zz[:, iy, ix], self.zc[:, iy, ix],
                        self.psi[it, ipsi, :, iy, ix]) for ix
                        in range(self.nx)] for iy in range(self.ny)] for ipsi
                        in range(self.npsi)] for it in range(self.nt)]
                psi = numpy.array(psi)
                self.psi[:, :, :, :, :] = numpy.moveaxis(psi, 4, 2)
                self.zc = self.zz

        else:

            # Model uses Cartesian height coordinate.
            self.zc = self.zz

    def write_data(self, path):
        """Write data into a NETCDF4 file."""

        rho = xarray.DataArray(self.psi[:, 0, :, :, :], dims = ("tt", "zz",
                "yy", "xx"), coords = {"tt": self.tt, "zz": self.zz[:, 0, 0],
                "yy": self.yy[0, :, 0], "xx": self.xx[0, 0, :]})
        uu = xarray.DataArray(self.psi[:, 1, :, :, :], dims = ("tt", "zz",
                "yy", "xx"), coords = {"tt": self.tt, "zz": self.zz[:, 0, 0],
                "yy": self.yy[0, :, 0], "xx": self.xx[0, 0, :]})
        vv = xarray.DataArray(self.psi[:, 2, :, :, :], dims = ("tt", "zz",
                "yy", "xx"), coords = {"tt": self.tt, "zz": self.zz[:, 0, 0],
                "yy": self.yy[0, :, 0], "xx": self.xx[0, 0, :]})
        ww = xarray.DataArray(self.psi[:, 3, :, :, :], dims = ("tt", "zz",
                "yy", "xx"), coords = {"tt": self.tt, "zz": self.zz[:, 0, 0],
                "yy": self.yy[0, :, 0], "xx": self.xx[0, 0, :]})
        pi = xarray.DataArray(self.psi[:, 4, :, :, :], dims = ("tt", "zz",
                "yy", "xx"), coords = {"tt": self.tt, "zz": self.zz[:, 0, 0],
                "yy": self.yy[0, :, 0], "xx": self.xx[0, 0, :]})
        zc = xarray.DataArray(self.zc, dims = ("zz", "yy", "xx"),
                coords = {"zz": self.zz[:, 0, 0], "yy": self.yy[0, :, 0],
                "xx": self.xx[0, 0, :]})
        self.data = xarray.Dataset({"rho": rho, "uu": uu, "vv": vv, "ww": ww,
                "pi": pi, "zc": zc})
        self.data.to_netcdf("".join((path, "data.nc")))
