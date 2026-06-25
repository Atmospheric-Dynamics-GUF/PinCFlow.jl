function add_contour_plot! end

function add_contour_plot!(figure::Figure, input::NamedTuple)
    (;
        colormap_name,
        columns,
        number,
        phi,
        row,
        significant_digits,
        title,
        x,
        x_label,
        y,
        y_label,
    ) = input

    sigdigits = significant_digits
    axis =
        Axis(figure[row, columns[1]]; title, xlabel = x_label, ylabel = y_label)
    (levels, colormap) =
        symmetric_contours(minimum(phi), maximum(phi); number, colormap_name)
    plot = contourf!(
        round.(x; sigdigits),
        round.(y; sigdigits),
        round.(phi; sigdigits);
        round.(levels; sigdigits),
        colormap,
    )
    tightlimits!(axis)
    Colorbar(
        figure[row, columns[2]],
        plot;
        ticks = round.(levels; sigdigits),
        tickformat = "{:8.1E}",
        label,
    )
    xlims!(minimum(x), maximum(x))
    ylims!(minimum(y), maximum(y))

    return
end
