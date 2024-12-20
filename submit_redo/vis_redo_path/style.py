import matplotlib.pyplot as pyplot
import matplotlib.colors as colors
import matplotlib.cm as cm
import cmasher
import numpy
import re

def symmetric_contours(minimum, maximum, number = 25, colormap = "default"):
  """Compute symmetric contour levels and crop the colormap accordingly."""

  # Compute contour levels.
  if minimum == - maximum or numpy.sign(minimum) == numpy.sign(maximum):
    levels = numpy.linspace(minimum, maximum, number)
  else:
    peak = max(abs(minimum), abs(maximum))
    factor = int(2.0 * peak / (maximum - minimum)) + 1
    if number % 2 and not factor % 2:
      factor += 1
    levels = numpy.linspace(- peak, peak, factor * number)
    if peak > - minimum:
      while levels[- number] > minimum:
        levels = numpy.linspace(- peak, peak, len(levels) - 2)
      levels = levels[(- number):]
    elif peak > maximum:
      while levels[number - 1] < maximum:
        levels = numpy.linspace(- peak, peak, len(levels) - 2)
      levels = levels[:number]

  # Get the colormap and peak.
  if isinstance(colormap, str):
    colormap = cm.get_cmap(colormap)
  peak = max(abs(levels[0]), abs(levels[- 1]))

  # Crop the colormap.
  if levels[0] > 0.0:
    colormap = colors.LinearSegmentedColormap.from_list("", \
        colormap(numpy.linspace(0.5, 1.0, 100)))
  elif levels[- 1] < 0.0:
    colormap = colors.LinearSegmentedColormap.from_list("", \
        colormap(numpy.linspace(0.0, 0.5, 100)))
  elif peak + levels[0] > 0.0:
    fraction = 0.5 * (peak + levels[0]) / peak
    colormap = colors.LinearSegmentedColormap.from_list("", \
        colormap(numpy.linspace(fraction, 1.0, 100)))
  elif peak - levels[- 1] > 0.0:
    fraction = 0.5 * (peak - levels[- 1]) / peak
    colormap = colors.LinearSegmentedColormap.from_list("", \
        colormap(numpy.linspace(0.0, 1.0 - fraction, 100)))

  # Return the results.
  return (levels, colormap)

# Define LaTeX preamble.
preamble = r"\usepackage{amsmath, amstext, amssymb, amsfonts, amsthm}" \
    r"\allowdisplaybreaks" \
    r"\usepackage[slantedGreek]{newtxmath}" \
    r"\renewcommand*\rmdefault{ptm}" \
    r"\renewcommand*\sfdefault{phv}" \
    r"\renewcommand*\ttdefault{lmtt}"

# Register default colormap.
bottom = cmasher.voltage(numpy.linspace(0.25, 1.0, 100))
top = cmasher.flamingo_r(numpy.linspace(0.0, 0.75, 100))
default = colors.LinearSegmentedColormap.from_list("default", \
    numpy.vstack((bottom, top)))
pyplot.register_cmap("default", default)

# Set plot parameters.
pyplot.style.use([entry for entry in pyplot.style.available if \
    re.search(r".*seaborn.*bright.*", entry)][0])
pyplot.rcParams["figure.constrained_layout.use"] = True
pyplot.rcParams["figure.figsize"] = (4.0, 3.0)
pyplot.rcParams["font.family"] = "serif"
pyplot.rcParams["image.cmap"] = "default"
pyplot.rcParams["legend.frameon"] = False
pyplot.rcParams["text.latex.preamble"] = preamble
pyplot.rcParams["text.usetex"] = True