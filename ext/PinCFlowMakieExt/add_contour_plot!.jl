"""
```julia
add_contour_plot!(figure::Figure, input::NamedTuple)
```

Add a contour plot to `figure` using the data in `input`.

# Arguments

  - `figure`: Figure to modify.

  - `input`: Data and specifications for the plot.

# See also

  - [`PinCFlow.symmetric_contours`](@ref)

  - [`PinCFlowMakieExt.tick_indices`](@ref)
"""
function add_contour_plot! end

@ivy function add_contour_plot!(figure::Figure, input::NamedTuple)
    (;
        background_color,
        colormap_name,
        color_tick_format,
        columns,
        label,
        number,
        phi,
        row,
        title,
        x,
        x_label,
        x_tick_format,
        y,
        y_label,
        y_tick_format,
    ) = input

    xmax = maximum(x)
    xmin = minimum(x)
    ymax = maximum(y)
    ymin = minimum(y)

    axis = Axis(
        figure[row, columns[1]];
        backgroundcolor = background_color,
        title,
        xlabel = x_label,
        xticks = xmin .+ [0.1, 0.5, 0.9] .* (xmax .- xmin),
        xtickformat = x_tick_format,
        ylabel = y_label,
        yticks = ymin .+ [0.1, 0.5, 0.9] .* (ymax .- ymin),
        ytickformat = y_tick_format,
    )
    (levels, colormap) =
        symmetric_contours(minimum(phi), maximum(phi); number, colormap_name)
    contourf!(x, y, phi; levels, colormap)
    tightlimits!(axis)
    Colorbar(
        figure[row, columns[2]];
        colormap = cgrad(colormap; categorical = true),
        label,
        limits = (levels[1], levels[end]),
        tickformat = color_tick_format,
        ticks = levels[tick_indices(levels)],
    )
    xlims!(xmin, xmax)
    ylims!(ymin, ymax)

    return
end
