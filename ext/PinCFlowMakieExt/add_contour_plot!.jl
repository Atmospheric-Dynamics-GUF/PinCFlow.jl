function add_contour_plot! end

function add_contour_plot!(figure::Figure, input::NamedTuple)
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

    axis = Axis(
        figure[row, columns[1]];
        backgroundcolor = background_color,
        title,
        xlabel = x_label,
        xtickformat = x_tick_format,
        ylabel = y_label,
        ytickformat = y_tick_format,
    )
    (levels, colormap) =
        symmetric_contours(minimum(phi), maximum(phi); number, colormap_name)
    plot = contourf!(x, y, phi; levels, colormap)
    tightlimits!(axis)
    Colorbar(
        figure[row, columns[2]],
        plot;
        ticks = levels,
        tickformat = color_tick_format,
        label,
    )
    xlims!(minimum(x), maximum(x))
    ylims!(minimum(y), maximum(y))

    return
end
