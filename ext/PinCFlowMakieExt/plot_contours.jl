function plot_contours(
    file::AbstractString,
    data::HDF5.File,
    variable::AbstractString,
    indices::Vararg{NTuple{4, <:Integer}};
    number::Integer = 10,
    colormap_name::Symbol = :seismic,
    label::AbstractString = "",
)
    set_visualization_theme!()

    # Set the grid.
    x = data["x"][:] ./ 1000
    y = data["y"][:] ./ 1000
    z = data["z"][:, :, :] ./ 1000
    (nx, ny, nz) = size(z)
    x = [xi for xi in x, j in 1:ny, k in 1:nz]
    y = [yj for i in 1:nx, yj in y, k in 1:nz]

    # Get the time.
    t = data["t"][:] ./ 3600

    # Get the variable.
    phi = data[variable][:, :, :, :]

    # Create the figure.
    figure = Figure()

    # Loop over outputs.
    row = 0
    for (i, j, k, n) in indices
        row += 1
        column = 0

        # Round the time.
        tn = round(t[n]; digits = 1)

        # Plot in the x-y plane.
        if nx > 1 && ny > 1
            column += 2
            @ivy zk = round(sum(z[:, :, k]) / length(z[:, :, k]); digits = 1)
            axis = Axis(
                figure[row, column - 1];
                title = L"t\approx%$tn\,\mathrm{h},\quad z\approx%$zk\,\mathrm{km}",
                xlabel = L"x\,[\mathrm{km}]",
                ylabel = L"y\,[\mathrm{km}]",
            )
            @ivy (levels, colormap) = symmetric_contours(
                minimum(phi[:, :, k, n]),
                maximum(phi[:, :, k, n]);
                number,
                colormap_name,
            )
            @ivy contours = contourf!(
                axis,
                x[:, :, k],
                y[:, :, k],
                phi[:, :, k, n];
                levels,
                colormap,
            )
            tightlimits!(axis)
            Colorbar(
                figure[row, column],
                contours;
                ticks = trunc.(levels; digits = 4),
                label,
            )
            xlims!(minimum(x), maximum(x))
            ylims!(minimum(y), maximum(y))
        end

        # Plot in the x-z plane.
        if nx > 1 && nz > 1
            column += 2
            @ivy yj = round(sum(y[:, j, :]) / length(y[:, j, :]); digits = 1)
            axis = Axis(
                figure[row, column - 1];
                title = L"t\approx%$tn\,\mathrm{h},\quad y\approx%$yj\,\mathrm{km}",
                xlabel = L"x\,[\mathrm{km}]",
                ylabel = L"z\,[\mathrm{km}]",
            )
            @ivy (levels, colormap) = symmetric_contours(
                minimum(phi[:, j, :, n]),
                maximum(phi[:, j, :, n]);
                number,
                colormap_name,
            )
            @ivy contours = contourf!(
                axis,
                x[:, j, :],
                z[:, j, :],
                phi[:, j, :, n];
                levels,
                colormap,
            )
            tightlimits!(axis)
            Colorbar(
                figure[row, column],
                contours;
                ticks = trunc.(levels; digits = 4),
                label,
            )
            xlims!(minimum(x), maximum(x))
            ylims!(minimum(z), maximum(z))
        end

        # Plot in the y-z plane.
        if ny > 1 && nz > 1
            column += 2
            @ivy xi = round(sum(x[i, :, :]) / length(x[i, :, :]); digits = 1)
            axis = Axis(
                figure[row, column - 1];
                title = L"t\approx%$tn\,\mathrm{h},\quad x\approx%$xi\,\mathrm{km}",
                xlabel = L"y\,[\mathrm{km}]",
                ylabel = L"z\,[\mathrm{km}]",
            )
            @ivy (levels, colormap) = symmetric_contours(
                minimum(phi[i, :, :, n]),
                maximum(phi[i, :, :, n]);
                number,
                colormap_name,
            )
            @ivy contours = contourf!(
                axis,
                y[i, :, :],
                z[i, :, :],
                phi[i, :, :, n];
                levels,
                colormap,
            )
            tightlimits!(axis)
            Colorbar(
                figure[row, column],
                contours;
                ticks = trunc.(levels; digits = 4),
                label,
            )
            xlims!(minimum(y), maximum(y))
            ylims!(minimum(z), maximum(z))
        end
    end

    # Resize, display and save the figure.
    resize_to_layout!(figure)
    display(figure)
    save(file, figure)

    return
end
