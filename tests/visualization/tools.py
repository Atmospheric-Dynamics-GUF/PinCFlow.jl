"""This Python module provides tools for working with PincFlow."""

import numpy
import numpy.linalg as linalg
import scipy.interpolate as interpolate
import matplotlib.pyplot as pyplot

class ModelOutput:
  """Import and transform model output data from PincFlow."""

  def __init__(self, path):
    """Import and format data."""

    # Read input file.
    with open("".join((path, "input.f90"))) as file:
      lines = file.readlines()

    # Remove white spaces, empty lines, comments and namelist titles.
    lines = [entry.split("!", 1)[0].strip() for entry in lines if \
        len(entry.strip()) > 1 and entry.strip()[0] not in ("!", "&")]

    # Get keys and values.
    keys = [entry.split("=", 1)[0].strip() for entry in lines]
    values = [entry.split("=", 1)[1].strip() for entry in lines]

    # Remove trailing commas.
    values = [entry.rstrip(",") if "," in entry else entry for entry in values]

    # Adjust boolean entries.
    values = [entry.replace(".true.", "True") if ".true." in entry else entry \
        for entry in values]
    values = [entry.replace(".false.", "False") if ".false." in entry else \
        entry for entry in values]

    # Evaluate entries.
    values = [eval(entry) if "=" not in entry else entry for entry in values]

    # Save to dictionary.
    parameters = {keys[nn]: values[nn] for nn in range(len(keys))}
    self.parameters = parameters

    # Set attributes.
    self.nx = parameters["sizeX"]
    self.ny = parameters["sizeY"]
    self.nz = parameters["sizeZ"]
    self.lx = parameters["lx_dim"][1] - parameters["lx_dim"][0]
    self.ly = parameters["ly_dim"][1] - parameters["ly_dim"][0]
    self.lz = parameters["lz_dim"][1] - parameters["lz_dim"][0]
    self.topography = parameters["topography"]
    self.topography_time = parameters["topographyTime"]
    self.output_type = parameters["outputType"]
    self.npsi = numpy.sum(parameters["varOut"])

    # Define grid spacing.
    self.dx = self.lx / self.nx
    self.dy = self.ly / self.ny
    self.dz = self.lz / self.nz

    # Define grid.
    self.xx = numpy.linspace(0.5 * self.dx, self.lx - 0.5 * self.dx, self.nx)
    self.yy = numpy.linspace(0.5 * self.dy, self.ly - 0.5 * self.dy, self.ny)
    self.zz = numpy.linspace(0.5 * self.dz, self.lz - 0.5 * self.dz, self.nz)
    (self.zz, self.yy, self.xx) = numpy.meshgrid(self.zz, self.yy, self.xx, \
        indexing = "ij")

    # Import data.
    psi = numpy.fromfile("".join((path, "pf_all.dat")), dtype = "float32")
    self.nt = int(len(psi) / self.nx / self.ny / self.nz / self.npsi)
    self.psi = numpy.reshape(psi, (self.nt, self.npsi, self.nz, self.ny, \
        self.nx))

    # Interpolate velocities to cell centers.
    if self.nx > 1:
      self.psi[:, 1, ..., 1:] = 0.5 * (self.psi[:, 1, ..., 1:] + self.psi[:, \
          1, ..., :(- 1)])
    if self.ny > 1:
      self.psi[:, 2, :, 1:] = 0.5 * (self.psi[:, 2, :, 1:] + self.psi[:, 2, :, \
          :(- 1)])
    if self.nz > 1:
      self.psi[:, 3, 1:] = 0.5 * (self.psi[:, 3, 1:] + self.psi[:, 3, :(- 1)])

    # Compute output times.
    self.tt = numpy.arange(self.nt, dtype = "float32")
    if self.output_type == "time":
      self.tt *= parameters["outputTimeDiff"]
    elif self.output_type == "timeStep":
      self.tt[1:] = numpy.loadtxt("".join((path, "dt.dat")), \
          dtype = "float32")[:, 1].cumsum()[(parameters["nOutput"] \
          - 1)::parameters["nOutput"]]

    # Define topography.
    if self.topography:
      hh = numpy.fromfile("".join((path, "topography.dat")), dtype = "float32")
      self.hh = numpy.reshape(hh, (self.ny, self.nx)) * numpy.ones((self.nt, \
          self.ny, self.nx))
      if self.topography_time > 0.0:
        self.hh[self.tt < self.topography_time] *= self.tt[self.tt \
            < self.topography_time, None, None] / self.topography_time

    # Import background fields.
    if self.parameters["model"] != "Boussinesq":
      pbar = numpy.fromfile("".join((path, "pStrat.dat")), dtype = "float32")
      thetabar = numpy.fromfile("".join((path, "thetaStrat.dat")), dtype \
          = "float32")
      rhobar = numpy.fromfile("".join((path, "rhoStrat.dat")), dtype \
          = "float32")
      n2bar = numpy.fromfile("".join((path, "bvsStrat.dat")), dtype \
          = "float32")
      if self.topography:
        self.pbar = pbar.reshape((self.nt, self.nz, self.ny, self.nx))
        self.thetabar = thetabar.reshape((self.nt, self.nz, self.ny, self.nx))
        self.rhobar = rhobar.reshape((self.nt, self.nz, self.ny, self.nx))
        self.n2bar = n2bar.reshape((self.nt, self.nz, self.ny, self.nx))
      else:
        self.pbar = pbar.reshape((self.nt, self.nz))
        self.thetabar = thetabar.reshape((self.nt, self.nz))
        self.rhobar = rhobar.reshape((self.nt, self.nz))
        self.n2bar = n2bar.reshape((self.nt, self.nz))

    # Import WKB data.
    if self.parameters["rayTracer"]:
      wkb = numpy.fromfile("".join((path, "pf_wkb_mean.dat")), dtype \
          = "float32")
      self.wkb = numpy.reshape(wkb, (self.nt, 13, self.nz, self.ny, self.nx))

  def transform(self, interpolation = False):
    """Transform and interpolate data."""

    if self.topography:

      # Define Jacobian.
      jj = (self.lz - self.hh) / self.lz

      # Define metric tensor elements.
      gg = numpy.zeros((2, self.nt, self.nz, self.ny, self.nx))
      gg[0, ..., 1:(- 1)] = (0.5 * (self.hh[:, None, :, 2:] - self.hh[:, None, \
          :, :(- 2)]) / self.dx * (self.zz[..., 1:(- 1)] - self.lz) / (self.lz \
          - self.hh[:, None, :, 1:(- 1)]))
      gg[1, :, :, 1:(- 1)] = (0.5 * (self.hh[:, None, 2:] - self.hh[:, None, \
          :(- 2)]) / self.dy * (self.zz[:, 1:(- 1)] - self.lz) / (self.lz \
          - self.hh[:, None, 1:(- 1)]))

      # Define metric tensor elements at the horizontal boundaries.
      if self.nx > 1:
        gg[0, ..., 0] = (0.5 * (self.hh[:, None, :, 1] - self.hh[:, None, :, \
            - 1]) / self.dx * (self.zz[..., 0] - self.lz) / (self.lz \
            - self.hh[:, None, :, 0]))
        gg[0, ..., - 1] = (0.5 * (self.hh[:, None, :, 0] - self.hh[:, None, :, \
            - 2]) / self.dx * (self.zz[..., - 1] - self.lz) / (self.lz \
            - self.hh[:, None, :, - 1]))
      if self.ny > 1:
        gg[1, :, :, 0] = (0.5 * (self.hh[:, None, 1] - self.hh[:, None, - 1]) \
            / self.dy * (self.zz[:, 0] - self.lz) / (self.lz - self.hh[:, \
            None, 0]))
        gg[1, :, :, - 1] = (0.5 * (self.hh[:, None, 0] - self.hh[:, None, \
            - 2]) / self.dy * (self.zz[:, - 1] - self.lz) / (self.lz \
            - self.hh[:, None, - 1]))

      # Compute Cartesian height.
      self.zc = (jj[:, None] * self.zz + self.hh[:, None]).copy()

      # Compute Cartesian vertical wind.
      self.psi[:, 3] = (jj[:, None] * self.psi[:, 3] - jj[:, None] * gg[0] \
          * self.psi[:, 1] - jj[:, None] * gg[1] * self.psi[:, 2])

      # Interpolate to Cartesian grid (NumPy is way faster than SciPy).
      if interpolation:
        for it in range(self.nt):
          for iy in range(self.ny):
            for ix in range(self.nx):
              for ipsi in range(self.npsi):
                psi = numpy.interp(self.zz[:, iy, ix], self.zc[it, :, iy, ix], \
                    self.psi[it, ipsi, :, iy, ix])
                self.psi[it, ipsi, :, iy, ix] = psi
        if self.parameters["model"] != "Boussinesq":
          for it in range(self.nt):
            for iy in range(self.ny):
              for ix in range(self.nx):
                pbar = numpy.interp(self.zz[:, iy, ix], self.zc[it, :, iy, \
                    ix], self.pbar[it, :, iy, ix])
                self.pbar[it, :, iy, ix] = pbar
                thetabar = numpy.interp(self.zz[:, iy, ix], self.zc[it, :, iy, \
                    ix], self.thetabar[it, :, iy, ix])
                self.thetabar[it, :, iy, ix] = thetabar
                rhobar = numpy.interp(self.zz[:, iy, ix], self.zc[it, :, iy, \
                    ix], self.rhobar[it, :, iy, ix])
                self.rhobar[it, :, iy, ix] = rhobar
                n2bar = numpy.interp(self.zz[:, iy, ix], self.zc[it, :, iy, \
                    ix], self.n2bar[it, :, iy, ix])
                self.n2bar[it, :, iy, ix] = n2bar
        if self.parameters["rayTracer"]:
          for it in range(self.nt):
            for iwkb in range(self.wkb.shape[1]):
              for iy in range(self.ny):
                for ix in range(self.nx):
                  wkb = numpy.interp(self.zz[:, iy, ix], self.zc[it, :, iy, \
                      ix], self.wkb[it, iwkb, :, iy, ix])
                  self.wkb[it, iwkb, :, iy, ix] = wkb
        self.zc = (numpy.ones((self.nt, self.nz, self.ny, self.nx)) \
            * self.zz).copy()

    else:

      # Model uses Cartesian height coordinate.
      self.zc = (numpy.ones((self.nt, self.nz, self.ny, self.nx)) \
          * self.zz).copy()

  def compute_wkb_topography(self, nxt, nyt, modes, regularization = 1.0):
    """Compute topographic input for WKB model."""

    # Define target grid.
    lxt = self.lx
    lyt = self.ly
    dxt = lxt / nxt
    dyt = lyt / nyt
    xt = numpy.linspace(0.5 * dxt, lxt - 0.5 * dxt, nxt)
    yt = numpy.linspace(0.5 * dyt, lyt - 0.5 * dyt, nyt)
    (yt, xt) = numpy.meshgrid(yt, xt, indexing = "ij")

    # Define interpolation grid.
    nxi = self.nx
    nyi = self.ny
    rx = max(2 * int(0.5 * nxi / nxt), 1)
    ry = max(2 * int(0.5 * nyi / nyt), 1)
    deltaix = int(0.5 * rx)
    deltaiy = int(0.5 * ry)
    dxi = dxt / rx
    dyi = dyt / ry
    lxi = nxi * dxi
    lyi = nyi * dyi
    xi = numpy.arange(xt[0, 0] - 0.5 * dxt, xt[0, 0] - 0.5 * dxt + lxi, dxi)
    yi = numpy.arange(yt[0, 0] - 0.5 * dyt, yt[0, 0] - 0.5 * dyt + lyi, dyi)
    (yi, xi) = numpy.meshgrid(yi, xi, indexing = "ij")

    # Find indices at target grid box centers.
    ixc = numpy.array([numpy.abs(xi[0] - xt[0, ix]).argmin() for ix in \
        range(nxt)])
    iyc = numpy.array([numpy.abs(yi[:, 0] - yt[iy, 0]).argmin() for iy in \
        range(nyt)])

    # Set grid box slices.
    if not deltaix:
      xslice = [slice(None) for ix in ixc]
    else:
      xslice = [slice(ix - deltaix, ix + deltaix) for ix in ixc]
    if not deltaiy:
      yslice = [slice(None) for iy in iyc]
    else:
      yslice = [slice(iy - deltaiy, iy + deltaiy) for iy in iyc]

    # Compute interpolated orography.
    if not deltaix:
      hi = interpolate.griddata((self.yy[0].flatten()), self.hh[- \
          1].flatten(), (yi.flatten()), fill_value = self.hh[- 1].mean(), \
          method = "cubic").reshape(nyi, nxi)
    elif not deltaiy:
      hi = interpolate.griddata((self.xx[0].flatten()), self.hh[- \
          1].flatten(), (xi.flatten()), fill_value = self.hh[- 1].mean(), \
          method = "cubic").reshape(nyi, nxi)
    else:
      hi = interpolate.griddata((self.yy[0].flatten(), self.xx[0].flatten()), \
          self.hh[- 1].flatten(), (yi.flatten(), xi.flatten()), fill_value \
          = self.hh[- 1].mean(), method = "cubic").reshape(nyi, nxi)

    # Compute mean orography on interpolation grid and subtract.
    hbar = hi.copy()
    aa = numpy.zeros((nyt, nxt))
    for iy in range(nyt):
      for ix in range(nxt):
        aa[iy, ix] = hi[yslice[iy], xslice[ix]].mean()
        hbar[yslice[iy], xslice[ix]] = aa[iy, ix]
    hi -= hbar

    # Compute spectral approximation for each grid box.
    hprime = numpy.zeros((nyi, nxi))
    bb = numpy.zeros((modes, nyt, nxt))
    kk = numpy.zeros((modes, nyt, nxt))
    ll = numpy.zeros((modes, nyt, nxt))
    pp = numpy.zeros((modes, nyt, nxt))
    kf = numpy.linspace(2.0 * numpy.pi / dxt, 2.0 * numpy.pi / dxi, rx)
    lf = numpy.linspace(2.0 * numpy.pi / dyt, 2.0 * numpy.pi / dyi, ry)
    pf = numpy.linspace(0.0, 0.5 * numpy.pi, rx * ry)
    (kf, lf, pf) = numpy.meshgrid(kf, lf, pf, indexing = "ij")
    kf = kf.flatten()
    lf = lf.flatten()
    pf = pf.flatten()
    for iy in range(nyt):
      for ix in range(nxt):
        yf = yi[yslice[iy], xslice[ix]].flatten()[:, None]
        xf = xi[yslice[iy], xslice[ix]].flatten()[:, None]
        hf = hi[yslice[iy], xslice[ix]].flatten()
        alpha = numpy.cos(kf * xf + lf * yf + pf)
        beta = linalg.inv(alpha.transpose() @ alpha + regularization \
            * numpy.eye(alpha.shape[1])) @ alpha.transpose() @ hf
        hprime[yslice[iy], xslice[ix]] = (alpha @ beta).reshape(ry, rx)
        maxima = numpy.argsort(numpy.abs(beta))[(- modes):]
        bb[:, iy, ix] = beta[maxima]
        kk[:, iy, ix] = kf[maxima]
        ll[:, iy, ix] = lf[maxima]
        pp[:, iy, ix] = pf[maxima]

    # Write input for the resolved topography.
    topography = numpy.array(aa, dtype = "float32").flatten()
    topography.tofile("topography.dat")

    # Write input for the unresolved topography.
    wkb_topography = numpy.array([numpy.abs(bb), numpy.abs(kk), \
        numpy.abs(ll)], dtype = "float32").flatten()
    wkb_topography.tofile("wkb_topography.dat")

    # Compare in a plot.
    peak = 2.0 * max(numpy.abs(self.hh).max(), numpy.abs(hprime).max())
    (figure, axes) = pyplot.subplots(1, 3, figsize = (8.0, 3.0))
    if nxt == 1:
      axes[0].plot(self.yy[0, :, 0], self.hh[- 1, :, 0])
      axes[1].plot(yi[:, 0], hbar[:, 0])
      axes[2].plot(yi[:, 0], hprime[:, 0] + hbar[:, 0])
      for index in range(len(axes)):
        axes[index].set_xlim(self.yy.min(), self.yy.max())
        axes[index].set_ylim(0.0, peak)
        axes[index].set_xlabel(r"$y \, \mathrm{\left[m\right]}$")
        axes[index].set_ylabel(r"$h \, \mathrm{\left[m\right]}$")
    elif nyt == 1:
      axes[0].plot(self.xx[0, 0], self.hh[- 1, 0])
      axes[1].plot(xi[0], hbar[0])
      axes[2].plot(xi[0], hprime[0] + hbar[0])
      for index in range(len(axes)):
        axes[index].set_xlim(self.xx.min(), self.xx.max())
        axes[index].set_ylim(0.0, peak)
        axes[index].set_xlabel(r"$x \, \mathrm{\left[m\right]}$")
        axes[index].set_ylabel(r"$h \, \mathrm{\left[m\right]}$")
    else:
      source = axes[0].pcolormesh(self.xx[0], self.yy[0], self.hh[- 1], vmin \
          = 0.0, vmax = peak, shading = "gouraud")
      mean = axes[1].pcolormesh(xi, yi, hbar, vmin = 0.0, vmax = peak, shading \
          = "gouraud")
      approximation = axes[2].pcolormesh(xi, yi, hprime + hbar, vmin = 0.0, \
          vmax = peak, shading = "gouraud")
      for index in range(len(axes)):
        axes[index].set_xlabel(r"$x \, \mathrm{\left[m\right]}$")
        axes[index].set_ylabel(r"$y \, \mathrm{\left[m\right]}$")
      figure.colorbar(source, ax = axes[0], label = r"$h \," \
          r"\mathrm{\left[m\right]}$")
      figure.colorbar(mean, ax = axes[1], label = r"$h \," \
          r"\mathrm{\left[m\right]}$")
      figure.colorbar(approximation, ax = axes[2], label = r"$h \," \
          r"\mathrm{\left[m\right]}$")
    axes[0].set_title("Source orography")
    axes[1].set_title("Mean orography")
    axes[2].set_title("Spectral approximation")
    figure.savefig("wkb_topography.pdf")
