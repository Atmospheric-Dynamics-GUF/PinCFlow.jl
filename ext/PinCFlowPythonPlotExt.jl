# Package extension for adding PythonPlot-based features to PinCFlow.jl.
module PinCFlowPythonPlotExt

using PythonPlot

"""
```julia
set_plot_style()
```

Configure PythonPlot.jl to use a preset plot style.
"""
function set_plot_style end

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

"""
```julia
symmetric_contours(
    minimum::AbstractFloat,
    maximum::AbstractFloat;
    number::Integer = 10,
    colormap_name::String = "seismic",
)::Tuple{<:LinRange{<:AbstractFloat, <:Integer}, <:Any}
```

Compute symmetric contours levels and return them and a correspondingly cropped colormap.

# Arguments

  - `minimum`: Smallest value to be plotted.

  - `maximum`: Largest value to be plotted.

# Keywords

  - `number`: Number of contour levels.

  - `colormap_name`: Name under which the chosen colormap is registered.
"""
function symmetric_contours end

function symmetric_contours(
    minimum::AbstractFloat,
    maximum::AbstractFloat;
    number::Integer = 10,
    colormap_name::String = "seismic",
)::Tuple{<:LinRange{<:AbstractFloat, <:Integer}, <:Any}

    # Compute contour levels.
    @ivy if minimum == -maximum ||
            sign(minimum) == sign(maximum) ||
            minimum == 0.0 ||
            maximum == 0.0
        levels = LinRange(minimum, maximum, number)
    else
        peak = max(abs(minimum), abs(maximum))
        factor = ceil(Int, 2.0 * peak / (maximum - minimum))
        if number % 2 > 0 && factor % 2 == 0
            factor += 1
        end
        levels = LinRange(-peak, peak, factor * number)
        if peak > -minimum
            while levels[end - number + 1] > minimum
                levels = LinRange(-peak, peak, length(levels) - 2)
            end
            levels = levels[(end - number + 1):end]
        elseif peak > maximum
            while levels[number] < maximum
                levels = LinRange(-peak, peak, length(levels) - 2)
            end
            levels = levels[1:number]
        end
    end

    # Get the colormap and peak.
    colormap = matplotlib.cm.get_cmap(colormap_name)
    @ivy peak = max(abs(levels[1]), abs(levels[end]))

    # Crop the colormap.
    @ivy if levels[1] > 0.0
        colormap = matplotlib.colors.LinearSegmentedColormap.from_list(
            "",
            colormap(LinRange(0.5, 1.0, 100)),
        )
    elseif levels[end] < 0.0
        colormap = matplotlib.colors.LinearSegmentedColormap.from_list(
            "",
            colormap(LinRange(0.0, 0.5, 100)),
        )
    elseif peak + levels[1] > 0.0
        fraction = 0.5 * (peak + levels[1]) / peak
        colormap = matplotlib.colors.LinearSegmentedColormap.from_list(
            "",
            colormap(LinRange(fraction, 1.0, 100)),
        )
    elseif peak - levels[end] > 0.0
        fraction = 0.5 * (peak - levels[end]) / peak
        colormap = matplotlib.colors.LinearSegmentedColormap.from_list(
            "",
            colormap(LinRange(0.0, 1.0 - fraction, 100)),
        )
    end

    return (levels, colormap)
end

end