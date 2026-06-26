function add_scatter_plot! end

function add_scatter_plot!(figure::Figure, input::NamedTuple)
    (;
        colormap_name,
        color_tick_format,
        columns,
        dx,
        dy,
        label,
        mask,
        number,
        phi,
        row,
        title,
        x,
        x_label,
        xmax,
        xmin,
        x_tick_format,
        y,
        y_label,
        ymax,
        ymin,
        y_tick_format,
    ) = input

    Axis(
        figure[row, columns[1]];
        title,
        xlabel = x_label,
        xtickformat = x_tick_format,
        ylabel = y_label,
        ytickformat = y_tick_format,
    )
    (levels, colormap) = symmetric_contours(
        minimum(phi[mask]),
        maximum(phi[mask]);
        number,
        colormap_name,
    )
    sorted = sortperm(abs.(phi[mask]))
    plot = scatter!(
        x[mask][sorted],
        y[mask][sorted];
        color = phi[mask][sorted],
        colormap = cgrad(colormap; categorical = true),
        marker = Rect,
        markersize = collect(zip(dx[mask][sorted], dy[mask][sorted])),
        markerspace = :data,
    )
    Colorbar(
        figure[row, columns[2]],
        plot;
        ticks = levels,
        tickformat = color_tick_format,
        label,
    )
    xlims!(xmin, xmax)
    ylims!(ymin, ymax)

    return
end
