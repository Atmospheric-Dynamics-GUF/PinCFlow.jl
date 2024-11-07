"""This Python module provides tools for working with PinCFlow."""

import re
import numpy
import numpy.fft as fft
import scipy.interpolate as interpolate
import matplotlib.pyplot as pyplot
import netCDF4
import style

class ModelOutput:
  """Import and transform model output data from PinCFlow."""

  def __init__(self, path):
    """Import and format data."""

    # Read namelists.
    with open("".join((path, "/namelists.txt"))) as file:
      lines = file.readlines()

    # Remove white spaces, empty lines, comments and namelist titles.
    lines = [entry.split("!", 1)[0].strip() for entry in lines if \
        len(entry.strip()) > 1 and entry.strip()[0] not in ("!", "&")]

    # Get keys and values.
    keys = [entry.split("=", 1)[0].strip().lower() for entry in lines]
    values = [entry.split("=", 1)[1].strip() for entry in lines]
    values = [entry.rstrip(",") if "," in entry else entry for entry in values]

    # Adjust boolean entries.
    values = [re.sub(r"\bT\b|\.true\.", "True", entry, flags = re.IGNORECASE) \
        if re.search(r"\bT\b|\.true\.", entry, flags = re.IGNORECASE) else \
        entry for entry in values]
    values = [re.sub(r"\bF\b|\.false\.", "False", entry, flags \
        = re.IGNORECASE) if re.search(r"\bF\b|\.false\.", entry, flags \
        = re.IGNORECASE) else entry for entry in values]

    # Adjust array entries.
    values = ["".join(("[", entry, "]")) if "," in entry or "*" in entry else \
        entry for entry in values]
    values = [re.sub(r"(\D)(\d+\*)([\d.eE+-]+)([^\d.eE+-]?)", r"\1\2[\3]\4", \
        entry) if re.search(r"(\D)(\d+\*)([\d.eE+-]+)([^\d.eE+-]*)", entry) \
        else entry for entry in values]

    # Evaluate entries.
    values = [eval(entry) if "=" not in entry else entry for entry in values]

    # Strip strings.
    values = [entry.strip() if isinstance(entry, str) else entry for entry in \
        values]

    # Flatten nested lists.
    values = [[[subentry] if not isinstance(subentry, list) else subentry for \
        subentry in entry] if isinstance(entry, list) else entry for entry in \
        values]
    values = [[subsubentry for subentry in entry for subsubentry in subentry] \
        if isinstance(entry, list) else entry for entry in values]

    # Concatenate array elements.
    values = [[values[nn] for nn in range(len(keys)) if keys[nn].split("(")[0] \
        == keys[mm].split("(")[0]] if "(" in keys[mm] else values[mm] for mm \
        in range(len(keys))]
    keys = [entry.split("(")[0] for entry in keys]

    # Save to dictionary.
    parameters = {keys[nn]: values[nn] for nn in range(len(keys))}
    self.parameters = parameters

    # Import data.
    data = netCDF4.Dataset("".join((path, "/pincflow_data_out.nc")))
    self.dimensions = data.dimensions
    self.variables = data.variables
    self.groups = data.groups

  def compute_wkb_topography(self, nxt, nyt):
    """Compute topographic input for WKB model."""

    # Get source grid and orography.
    lxs = self.parameters["lx_dim"][1] - self.parameters["lx_dim"][0]
    lys = self.parameters["ly_dim"][1] - self.parameters["ly_dim"][0]
    nxs = self.dimensions["x"].size
    nys = self.dimensions["y"].size
    dxs = lxs / nxs
    dys = lys / nys
    (ys, xs) = numpy.meshgrid(self.variables["y"], self.variables["x"], \
        indexing = "ij")
    hs = self.variables["h"][- 1]

    # Define target grid.
    lxt = lxs
    lyt = lys
    dxt = lxt / nxt
    dyt = lyt / nyt
    xt = numpy.linspace(0.5 * dxt, lxt - 0.5 * dxt, nxt)
    yt = numpy.linspace(0.5 * dyt, lyt - 0.5 * dyt, nyt)
    (yt, xt) = numpy.meshgrid(yt, xt, indexing = "ij")

    # Define interpolation grid.
    lxi = lxt
    lyi = lyt
    rx = max(nxs // nxt, 1)
    ry = max(nys // nyt, 1)
    dxi = dxt / rx
    dyi = dyt / ry
    nxi = nxt * rx
    nyi = nyt * ry
    xi = numpy.linspace(0.5 * dxi, lxi - 0.5 * dxi, nxi)
    yi = numpy.linspace(0.5 * dyi, lyi - 0.5 * dyi, nyi)
    (yi, xi) = numpy.meshgrid(yi, xi, indexing = "ij")

    # Set grid box slices.
    xslice = [slice(rx * ix, rx * (ix + 1)) for ix in range(nxt)]
    yslice = [slice(ry * iy, ry * (iy + 1)) for iy in range(nyt)]

    # Compute interpolated orography.
    if nxt == 1 and nyt == 1:
      hi = hs.copy()
    elif nxt == 1:
      hi = interpolate.griddata((ys.flatten()), hs.flatten(), (yi.flatten()), \
          method = "cubic").reshape(nyi, nxi)
    elif nyt == 1:
      hi = interpolate.griddata((xs.flatten()), hs.flatten(), (xi.flatten()), \
          method = "cubic").reshape(nyi, nxi)
    else:
      hi = interpolate.griddata((ys.flatten(), xs.flatten()), hs.flatten(), \
          (yi.flatten(), xi.flatten()), method = "cubic").reshape(nyi, nxi)

    # Compute mean orography on interpolation grid and subtract.
    hbar = numpy.zeros((nyt, nxt))
    hprime = hi.copy()
    for iy in range(nyt):
      for ix in range(nxt):
        hbar[iy, ix] = hi[yslice[iy], xslice[ix]].mean()
        hprime[yslice[iy], xslice[ix]] -= hbar[iy, ix]

    # Compute RFFT for each grid box (the factor two is necessary because
    # the output is truncated, the other factor takes care of the
    # normalization).
    ll = fft.fftfreq(ry)[:, None, None, None]
    kk = fft.rfftfreq(rx)[None, :, None, None]
    ll = ll * numpy.ones((1, kk.shape[1], nyt, nxt))
    kk = kk * numpy.ones((ll.shape[0], 1, nyt, nxt))
    htilde = numpy.zeros((ll.shape[0], kk.shape[1], nyt, nxt), dtype \
        = "complex64")
    for iy in range(nyt):
      for ix in range(nxt):
        htilde[..., iy, ix] = 2.0 * fft.rfft2(hprime[yslice[iy], xslice[ix]]) \
            / ry / rx

    # Flatten and remove zero-mode.
    htilde = htilde.flatten()[1:]
    kk = kk.flatten()[1:]
    ll = ll.flatten()[1:]

    # Write input for the resolved topography.
    topography = numpy.array(hbar, dtype = "float32").flatten()
    topography.tofile("topography.dat")

    # Write input for the unresolved topography.
    wkb_topography = numpy.array([htilde, kk, ll], dtype \
        = "complex64").flatten()
    wkb_topography.tofile("wkb_topography.dat")

    # Compare in a plot.
    (figure, axes) = pyplot.subplots(2, 2, figsize = (8.0, 6.0))
    if nxt == 1:
      axes[0, 0].plot(ys[:, 0], hs[:, 0])
      axes[0, 1].plot(yi[:, 0], hi[:, 0])
      axes[1, 0].plot(yt[:, 0], hbar[:, 0])
      axes[1, 1].plot(yi[:, 0], hprime[:, 0])
      axes[0, 0].set_ylabel(r"$h \, \mathrm{\left[m\right]}$")
      axes[0, 1].set_ylabel(r"$h_\mathrm{i} \, \mathrm{\left[m\right]}$")
      axes[1, 0].set_ylabel(r"$h_\mathrm{m} \, \mathrm{\left[m\right]}$")
      axes[1, 1].set_ylabel(r"$h_\mathrm{i} - h_\mathrm{m} \," \
          r"\mathrm{\left[m\right]}$")
      for entry in axes.flatten():
        entry.set_xlabel(r"$y \, \mathrm{\left[m\right]}$")
    elif nyt == 1:
      axes[0, 0].plot(xs[0], hs[0])
      axes[0, 1].plot(xi[0], hi[0])
      axes[1, 0].plot(xt[0], hbar[0])
      axes[1, 1].plot(xi[0], hprime[0])
      axes[0, 0].set_ylabel(r"$h \, \mathrm{\left[m\right]}$")
      axes[0, 1].set_ylabel(r"$h_\mathrm{i} \, \mathrm{\left[m\right]}$")
      axes[1, 0].set_ylabel(r"$h_\mathrm{m} \, \mathrm{\left[m\right]}$")
      axes[1, 1].set_ylabel(r"$h_\mathrm{i} - h_\mathrm{m} \," \
          r"\mathrm{\left[m\right]}$")
      for entry in axes.flatten():
        entry.set_xlabel(r"$x \, \mathrm{\left[m\right]}$")
    else:
      minimum = min(hs.min(), hi.min(), hbar.min(), hprime.min())
      maximum = max(hs.max(), hi.max(), hbar.max(), hprime.max())
      (levels, colormap) = style.symmetric_contours(minimum, maximum)
      source = axes[0, 0].contourf(xs, ys, hs, levels, cmap = colormap)
      interpolation = axes[0, 1].contourf(xi, yi, hi, levels, cmap = colormap)
      mean = axes[1, 0].contourf(xt, yt, hbar, levels, cmap = colormap)
      deviation = axes[1, 1].contourf(xi, yi, hprime, levels, cmap = colormap)
      figure.colorbar(source, ax = axes[0, 0], label = r"$h \," \
          r"\mathrm{\left[m\right]}$")
      figure.colorbar(interpolation, ax = axes[0, 1], label = r"$h_\mathrm{i}" \
          r"\, \mathrm{\left[m\right]" r"}$")
      figure.colorbar(mean, ax = axes[1, 0], label = r"$h_\mathrm{m} \," \
          r"\mathrm{\left[m\right]}$")
      figure.colorbar(deviation, ax = axes[1, 1], label = r"$h_\mathrm{i} -" \
          r"h_\mathrm{m} \, \mathrm{\left[m\right]}$")
      for entry in axes.flatten():
        entry.set_xlabel(r"$x \, \mathrm{\left[m\right]}$")
        entry.set_ylabel(r"$y \, \mathrm{\left[m\right]}$")
    axes[0, 0].set_title("Source orography")
    axes[0, 1].set_title("Interpolated orography")
    axes[1, 0].set_title("Mean orography")
    axes[1, 1].set_title("Deviation from mean")
    figure.savefig("wkb_topography.pdf", transparent = True)
    figure.savefig("wkb_topography.png")