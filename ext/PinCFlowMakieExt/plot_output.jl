function plot_output(
    file::AbstractString,
    data::HDF5.File,
    fields::Vararg{
        Tuple{<:AbstractString, <:Integer, <:Integer, <:Integer, <:Integer},
    };
    number::Integer = 10,
    colormap_name::Symbol = :seismic,
    space_unit::AbstractString = "km",
    time_unit::AbstractString = "h",
)
    set_visualization_theme!()

    # Store the ray-volume property names.
    ray_volume_properties = (
        "xr",
        "yr",
        "zr",
        "dxr",
        "dyr",
        "dzr",
        "kr",
        "lr",
        "mr",
        "dkr",
        "dlr",
        "dmr",
        "nr",
    )

    # Set the space unit factor.
    if space_unit == "km"
        space_unit_factor = 1000
    elseif space_unit == "m"
        space_unit_factor = 1
    else
        error("Error: Unknown space unit!")
    end

    # Set the time unit factor.
    if time_unit == "d"
        time_unit_factor = 86400
    elseif time_unit == "h"
        time_unit_factor = 3600
    elseif time_unit == "min"
        time_unit_factor = 60
    elseif time_unit == "s"
        time_unit_factor = 1
    else
        error("Error: Unknown time unit!")
    end

    # Set the grid.
    x = data["x"][:] ./ space_unit_factor
    y = data["y"][:] ./ space_unit_factor
    z = data["z"][:, :, :] ./ space_unit_factor
    (nx, ny, nz) = size(z)
    x = [xi for xi in x, j in 1:ny, k in 1:nz]
    y = [yj for i in 1:nx, yj in y, k in 1:nz]

    # Get the time.
    t = data["t"][:] ./ time_unit_factor

    # Create the figure.
    figure = Figure()

    # Loop over outputs.
    row = 0
    for (variable, i, j, k, n) in fields
        row += 1
        column = 0

        # Round the time.
        tn = round(t[n]; digits = 1)

        if variable in ("wavespectrum", )
            column += 2

             kp = data["kp"][:] .* space_unit_factor
             m = data["m"][:] .* space_unit_factor
             phi = data[variable][:, :, :, :, :, n]
             nkp = length(kp)
             nm = length(m)
             kp = [kpi for kpi in kp, mi in 1:nm]
             m = [mi for kpi in 1:nkp, mi in m]
             # Get the label.
            label = LaTeXString(attrs(data[variable])["label"])
            # Plot in the kp-m plane.
            if nkp > 1 && nm > 1
                column += 2
                @ivy zk =
                    round(sum(z[:, :, k]) / length(z[:, :, k]); digits = 1)
                @ivy yj =
                    round(sum(y[:, j, :]) / length(y[:, j, :]); digits = 1)
                @ivy xi =
                    round(sum(x[i, :, :]) / length(x[i, :, :]); digits = 1)
                axis = Axis(
                    figure[row, column - 1];
                    title = L"t\approx%$tn\ \mathrm{%$time_unit}, \quad x\approx%$xi\ \mathrm{%$space_unit}, 
                     \quad y\approx%$yj\ \mathrm{%$space_unit}, \quad z\approx%$zk\ \mathrm{%$space_unit}",
                    xlabel = L"kp \ [\mathrm{%$space_unit}^{-1}]",
                    ylabel = L"m \ [\mathrm{%$space_unit}^{-1}]" ,
                )
                @ivy (levels, colormap) = symmetric_contours(
                    minimum(phi[i, j, k, :, :]),
                    maximum(phi[i, j, k, :, :]);
                    number,
                    colormap_name,
                )
                @ivy plot = contourf!(
                    kp,
                    m,
                    phi[i,j,k, :, :];
                    levels,
                    colormap,
                )
                tightlimits!(axis)
                Colorbar(
                    figure[row, column],
                    plot;
                    ticks = levels,
                    tickformat = "{:9.2E}",
                    label,
                    
                )
                xlims!(minimum(kp), maximum(kp))
                ylims!(minimum(m), maximum(m))
            end

        elseif variable in ray_volume_properties
            # Get the ray-volume data.
            xr = data["xr"][:, :, :, :, n] ./ space_unit_factor
            yr = data["yr"][:, :, :, :, n] ./ space_unit_factor
            zr = data["zr"][:, :, :, :, n] ./ space_unit_factor
            dxr = data["dxr"][:, :, :, :, n] ./ space_unit_factor
            dyr = data["dyr"][:, :, :, :, n] ./ space_unit_factor
            dzr = data["dzr"][:, :, :, :, n] ./ space_unit_factor
            nr = data["nr"][:, :, :, :, n]
            phi = data[variable][:, :, :, :, n]

            # Get the label.
            label = LaTeXString(attrs(data[variable])["label"])

            # Plot in the x-y plane.
            if nx > 1 && ny > 1
                column += 2
                @ivy zk =
                    round(sum(z[:, :, k]) / length(z[:, :, k]); digits = 1)
                Axis(
                    figure[row, column - 1];
                    title = L"t\approx%$tn\ \mathrm{%$time_unit},\quad z\approx%$zk\ \mathrm{%$space_unit}",
                    xlabel = L"x_r\ [\mathrm{%$space_unit}]",
                    ylabel = L"y_r\ [\mathrm{%$space_unit}]",
                )
                @ivy nonzero = nr[:, :, :, k] .!= 0
                @ivy (levels, colormap) = symmetric_contours(
                    minimum(phi[:, :, :, k][nonzero]),
                    maximum(phi[:, :, :, k][nonzero]);
                    number,
                    colormap_name,
                )
                @ivy plot = scatter!(
                    xr[:, :, :, k][nonzero],
                    yr[:, :, :, k][nonzero];
                    color = phi[:, :, :, k][nonzero],
                    colormap = cgrad(colormap; categorical = true),
                    marker = Rect,
                    markersize = collect(
                        zip(dxr[:, :, :, k][nonzero], dyr[:, :, :, k][nonzero]),
                    ),
                    markerspace = :data,
                )
                Colorbar(
                    figure[row, column],
                    plot;
                    ticks = levels,
                    tickformat = "{:9.2E}",
                    label,
                )
                xlims!(minimum(x), maximum(x))
                ylims!(minimum(y), maximum(y))
            end

            # Plot in the x-z plane.
            if nx > 1 && nz > 1
                column += 2
                @ivy yj =
                    round(sum(y[:, j, :]) / length(y[:, j, :]); digits = 1)
                Axis(
                    figure[row, column - 1];
                    title = L"t\approx%$tn\ \mathrm{%$time_unit},\quad y\approx%$yj\ \mathrm{%$space_unit}",
                    xlabel = L"x_r\ [\mathrm{%$space_unit}]",
                    ylabel = L"z_r\ [\mathrm{%$space_unit}]",
                )
                @ivy nonzero = phi[:, :, j, :] .!= 0
                @ivy (levels, colormap) = symmetric_contours(
                    minimum(phi[:, :, j, :][nonzero]),
                    maximum(phi[:, :, j, :][nonzero]);
                    number,
                    colormap_name,
                )
                @ivy plot = scatter!(
                    xr[:, :, j, :][nonzero],
                    zr[:, :, j, :][nonzero];
                    color = phi[:, :, j, :][nonzero],
                    colormap = cgrad(colormap; categorical = true),
                    marker = Rect,
                    markersize = collect(
                        zip(dxr[:, :, j, :][nonzero], dzr[:, :, j, :][nonzero]),
                    ),
                    markerspace = :data,
                )
                Colorbar(
                    figure[row, column],
                    plot;
                    ticks = levels,
                    tickformat = "{:9.2E}",
                    label,
                )
                xlims!(minimum(x), maximum(x))
                ylims!(minimum(z), maximum(z))
            end

            # Plot in the y-z plane.
            if ny > 1 && nz > 1
                column += 2
                @ivy xi =
                    round(sum(x[i, :, :]) / length(x[i, :, :]); digits = 1)
                Axis(
                    figure[row, column - 1];
                    title = L"t\approx%$tn\ \mathrm{%$time_unit},\quad x\approx%$xi\ \mathrm{%$space_unit}",
                    xlabel = L"y_r\ [\mathrm{%$space_unit}]",
                    ylabel = L"z_r\ [\mathrm{%$space_unit}]",
                )
                @ivy nonzero = phi[:, i, :, :] .!= 0
                @ivy (levels, colormap) = symmetric_contours(
                    minimum(phi[:, i, :, :][nonzero]),
                    maximum(phi[:, i, :, :][nonzero]);
                    number,
                    colormap_name,
                )
                @ivy plot = scatter!(
                    yr[:, i, :, :][nonzero],
                    zr[:, i, :, :][nonzero];
                    color = phi[:, i, :, :][nonzero],
                    colormap = cgrad(colormap; categorical = true),
                    marker = Rect,
                    markersize = collect(
                        zip(dyr[:, i, :, :][nonzero], dzr[:, i, :, :][nonzero]),
                    ),
                    markerspace = :data,
                )
                Colorbar(
                    figure[row, column],
                    plot;
                    ticks = levels,
                    tickformat = "{:9.2E}",
                    label,
                )
                xlims!(minimum(y), maximum(y))
                ylims!(minimum(z), maximum(z))
            end
        else
            # Get the variable.
            phi = data[variable][:, :, :, n]

            # Get the label.
            label = LaTeXString(attrs(data[variable])["label"])

            # Plot in the x-y plane.
            if nx > 1 && ny > 1
                column += 2
                @ivy zk =
                    round(sum(z[:, :, k]) / length(z[:, :, k]); digits = 1)
                axis = Axis(
                    figure[row, column - 1];
                    title = L"t\approx%$tn\ \mathrm{%$time_unit},\quad z\approx%$zk\ \mathrm{%$space_unit}",
                    xlabel = L"x\ [\mathrm{%$space_unit}]",
                    ylabel = L"y\ [\mathrm{%$space_unit}]",
                )
                @ivy (levels, colormap) = symmetric_contours(
                    minimum(phi[:, :, k]),
                    maximum(phi[:, :, k]);
                    number,
                    colormap_name,
                )
                @ivy plot = contourf!(
                    x[:, :, k],
                    y[:, :, k],
                    phi[:, :, k];
                    levels,
                    colormap,
                )
                tightlimits!(axis)
                Colorbar(
                    figure[row, column],
                    plot;
                    ticks = levels,
                    tickformat = "{:9.2E}",
                    label,
                )
                xlims!(minimum(x), maximum(x))
                ylims!(minimum(y), maximum(y))
            end

            # Plot in the x-z plane.
            if nx > 1 && nz > 1
                column += 2
                @ivy yj =
                    round(sum(y[:, j, :]) / length(y[:, j, :]); digits = 1)
                axis = Axis(
                    figure[row, column - 1];
                    backgroundcolor = :black,
                    title = L"t\approx%$tn\ \mathrm{%$time_unit},\quad y\approx%$yj\ \mathrm{%$space_unit}",
                    xlabel = L"x\ [\mathrm{%$space_unit}]",
                    ylabel = L"z\ [\mathrm{%$space_unit}]",
                )
                @ivy (levels, colormap) = symmetric_contours(
                    minimum(phi[:, j, :]),
                    maximum(phi[:, j, :]);
                    number,
                    colormap_name,
                )
                @ivy plot = contourf!(
                    x[:, j, :],
                    z[:, j, :],
                    phi[:, j, :];
                    levels,
                    colormap,
                )
                tightlimits!(axis)
                Colorbar(
                    figure[row, column],
                    plot;
                    ticks = levels,
                    tickformat = "{:9.2E}",
                    label,
                )
                xlims!(minimum(x), maximum(x))
                ylims!(minimum(z), maximum(z))
            end

            # Plot in the y-z plane.
            if ny > 1 && nz > 1
                column += 2
                @ivy xi =
                    round(sum(x[i, :, :]) / length(x[i, :, :]); digits = 1)
                axis = Axis(
                    figure[row, column - 1];
                    backgroundcolor = :black,
                    title = L"t\approx%$tn\ \mathrm{%$time_unit},\quad x\approx%$xi\ \mathrm{%$space_unit}",
                    xlabel = L"y\ [\mathrm{%$space_unit}]",
                    ylabel = L"z\ [\mathrm{%$space_unit}]",
                )
                @ivy (levels, colormap) = symmetric_contours(
                    minimum(phi[i, :, :]),
                    maximum(phi[i, :, :]);
                    number,
                    colormap_name,
                )
                @ivy plot = contourf!(
                    y[i, :, :],
                    z[i, :, :],
                    phi[i, :, :];
                    levels,
                    colormap,
                )
                tightlimits!(axis)
                Colorbar(
                    figure[row, column],
                    plot;
                    ticks = levels,
                    tickformat = "{:9.2E}",
                    label,
                )
                xlims!(minimum(y), maximum(y))
                ylims!(minimum(z), maximum(z))
            end
        end
    end

    # Resize, display and save the figure.
    resize_to_layout!(figure)
    display(figure)
    save(file, figure)

    return
end
