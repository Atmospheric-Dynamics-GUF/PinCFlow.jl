import PinCFlow: set_plot_style
"""
```julia
set_plot_style()
```

Configure PythonPlot.jl to use a preset plot style.
"""
function set_plot_style()
    matplotlib.style.use(
        [
            entry for entry in matplotlib.style.available if
            occursin(r".*seaborn.*bright.*", string(entry))
        ][1],
    )
    matplotlib.rcParams["figure.autolayout"] = true
    matplotlib.rcParams["figure.figsize"] = (4.0, 3.0)
    matplotlib.rcParams["figure.dpi"] = 500
    matplotlib.rcParams["font.family"] = "serif"
    matplotlib.rcParams["image.cmap"] = "seismic"
    matplotlib.rcParams["legend.frameon"] = false
    matplotlib.rcParams["text.usetex"] = true
    matplotlib.rcParams["text.latex.preamble"] =
        "\\usepackage{amsmath, amstext, amssymb, amsfonts, amsthm}" *
        "\\allowdisplaybreaks" *
        # "\\usepackage[slantedGreek]{newtxmath}" *
        "\\renewcommand*\\rmdefault{ptm}" *
        "\\renewcommand*\\sfdefault{phv}" *
        "\\renewcommand*\\ttdefault{lmtt}"
    return
end