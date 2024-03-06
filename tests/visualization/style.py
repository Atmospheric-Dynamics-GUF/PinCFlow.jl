import matplotlib.pyplot as pyplot

# Define LaTeX preamble.
preamble = (r"\usepackage{amsmath, amstext, amssymb, amsfonts, amsthm}" \
    r"\allowdisplaybreaks" \
    r"\let\Gamma\varGamma" \
    r"\let\Delta\varDelta" \
    r"\let\Theta\varTheta" \
    r"\let\Lambda\varLambda" \
    r"\let\Xi\varXi" \
    r"\let\Pi\varPi" \
    r"\let\Sigma\varSigma" \
    r"\let\Upsilon\varUpsilon" \
    r"\let\Phi\varPhi" \
    r"\let\Psi\varPsi" \
    r"\let\Omega\varOmega" \
    r"\renewcommand*\rmdefault{ptm}" \
    r"\renewcommand*\sfdefault{phv}" \
    r"\renewcommand*\ttdefault{lmtt}")

# Set plot parameters.
pyplot.style.use("seaborn-v0_8-bright")
pyplot.rcParams["figure.constrained_layout.use"] = True
pyplot.rcParams["figure.figsize"] = (4.0, 3.0)
pyplot.rcParams["font.family"] = "serif"
pyplot.rcParams["image.cmap"] = "plasma"
pyplot.rcParams["legend.frameon"] = False
pyplot.rcParams["text.latex.preamble"] = preamble
pyplot.rcParams["text.usetex"] = True
