import matplotlib.pyplot as pyplot

# Define LaTeX preamble.
preamble = (r"\usepackage{amsmath, amstext, amssymb, amsfonts, amsthm}"
        r"\allowdisplaybreaks"
        r"\usepackage{bm}"
        r"\makeatletter"
        r"\let\alloc@latex\alloc@"
        r"\def\alloc@#1#2#3#4{\newcount}"
        r"\makeatother"
        r"\usepackage[lite, slantedGreek, compatiblegreek, "
        r"straightbraces]{mtpro2}"
        r"\makeatletter"
        r"\let\alloc@\alloc@latex"
        r"\makeatother"
        r"\DeclareSymbolFont{letters}{OML}{cmm}{m}{it}"
        r"\DeclareMathSymbol{0}{\mathalpha}{operators}{`0}"
        r"\DeclareMathSymbol{1}{\mathalpha}{operators}{`1}"
        r"\DeclareMathSymbol{2}{\mathalpha}{operators}{`2}"
        r"\DeclareMathSymbol{3}{\mathalpha}{operators}{`3}"
        r"\DeclareMathSymbol{4}{\mathalpha}{operators}{`4}"
        r"\DeclareMathSymbol{5}{\mathalpha}{operators}{`5}"
        r"\DeclareMathSymbol{6}{\mathalpha}{operators}{`6}"
        r"\DeclareMathSymbol{7}{\mathalpha}{operators}{`7}"
        r"\DeclareMathSymbol{8}{\mathalpha}{operators}{`8}"
        r"\DeclareMathSymbol{9}{\mathalpha}{operators}{`9}"
        r"\renewcommand*\rmdefault{ptm}"
        r"\renewcommand*\sfdefault{phv}"
        r"\renewcommand*\ttdefault{lmtt}")

# Set plot parameters.
pyplot.style.use("seaborn-bright")
pyplot.rcParams["figure.constrained_layout.use"] = True
pyplot.rcParams["figure.figsize"] = (4.0, 3.0)
pyplot.rcParams["font.family"] = "serif"
pyplot.rcParams["image.cmap"] = "plasma"
pyplot.rcParams["legend.frameon"] = False
pyplot.rcParams["text.latex.preamble"] = preamble
pyplot.rcParams["text.usetex"] = True