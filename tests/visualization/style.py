import matplotlib.pyplot as pyplot

# Define LaTeX preamble.
preamble = r"\usepackage{amsmath, amstext, amssymb, amsfonts, amsthm}" \
    r"\allowdisplaybreaks" \
    r"\usepackage[slantedGreek]{newtxmath}" \
    r"\renewcommand*\rmdefault{ptm}" \
    r"\renewcommand*\sfdefault{phv}" \
    r"\renewcommand*\ttdefault{lmtt}"

# Set plot parameters.
pyplot.style.use("seaborn-v0_8-bright")
pyplot.rcParams["figure.constrained_layout.use"] = True
pyplot.rcParams["figure.figsize"] = (4.0, 3.0)
pyplot.rcParams["font.family"] = "serif"
pyplot.rcParams["image.cmap"] = "plasma"
pyplot.rcParams["legend.frameon"] = False
pyplot.rcParams["text.latex.preamble"] = preamble
pyplot.rcParams["text.usetex"] = True
