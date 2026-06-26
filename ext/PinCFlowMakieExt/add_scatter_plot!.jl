function add_scatter_plot! end

@ivy function add_scatter_plot!(figure::Figure, input::NamedTuple)
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
        xticks = xmin .+ [0.1, 0.5, 0.9] .* (xmax .- xmin),
        xtickformat = x_tick_format,
        ylabel = y_label,
        yticks = ymin .+ [0.1, 0.5, 0.9] .* (ymax .- ymin),
        ytickformat = y_tick_format,
    )
    (levels, colormap) = symmetric_contours(
        minimum(phi[mask]; init = 0.0),
        maximum(phi[mask]; init = 0.0);
        number,
        colormap_name,
    )
    sorted = sortperm(abs.(phi[mask]))
    scatter!(
        x[mask][sorted],
        y[mask][sorted];
        color = phi[mask][sorted],
        colormap = cgrad(colormap; categorical = true),
        marker = Rect,
        markersize = collect(zip(dx[mask][sorted], dy[mask][sorted])),
        markerspace = :data,
    )
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
