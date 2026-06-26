function add_scatter_plot! end

function add_scatter_plot!(figure::Figure, input::NamedTuple)
    (;
        colormap_name,
        columns,
        dx,
        dy,
        label,
        mask,
        number,
        phi,
        row,
        significant_digits,
        title,
        x,
        x_label,
        xmax,
        xmin,
        y,
        y_label,
        ymax,
        ymin,
    ) = input

    sigdigits = significant_digits
    Axis(figure[row, columns[1]]; title, xlabel = x_label, ylabel = y_label)
    (levels, colormap) = symmetric_contours(
        minimum(phi[mask]),
        maximum(phi[mask]);
        number,
        colormap_name,
    )
    plot = scatter!(
        round.(x[mask]; sigdigits),
        round.(y[mask]; sigdigits);
        color = round.(phi[mask]; sigdigits),
        colormap = cgrad(colormap; categorical = true),
        marker = Rect,
        markersize = collect(
            zip(round.(dx[mask]; sigdigits), round.(dy[mask]; sigdigits)),
        ),
        markerspace = :data,
    )
    Colorbar(
        figure[row, columns[2]],
        plot;
        ticks = round.(levels; sigdigits),
        tickformat = "{:.$(sigdigits - 1)E}",
        label,
    )
    xlims!(xmin, xmax)
    ylims!(ymin, ymax)

    return
end
