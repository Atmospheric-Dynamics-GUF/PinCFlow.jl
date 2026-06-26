@ivy function plot_output(
    file::AbstractString,
    data::HDF5.File,
    fields::Vararg{
        Tuple{<:AbstractString, <:Integer, <:Integer, <:Integer, <:Integer},
    };
    colormap_name::Symbol = :seismic,
    significant_digits::Integer = 3,
    number::Integer = 10,
    space_unit::Symbol = :km,
    time_unit::Symbol = :h,
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
    if space_unit === :km
        space_unit_factor = 1000
    elseif space_unit === :m
        space_unit_factor = 1
    else
        error("Error: Unknown space unit!")
    end

    # Set the time unit factor.
    if time_unit === :d
        time_unit_factor = 86400
    elseif time_unit === :h
        time_unit_factor = 3600
    elseif time_unit === :min
        time_unit_factor = 60
    elseif time_unit === :s
        time_unit_factor = 1
    else
        error("Error: Unknown time unit!")
    end

    # Set the significant digits.
    sigdigits = significant_digits

    # Set the grid.
    x = round.(data["x"][:]; sigdigits) ./ space_unit_factor
    y = round.(data["y"][:]; sigdigits) ./ space_unit_factor
    z = round.(data["z"][:, :, :]; sigdigits) ./ space_unit_factor
    (nx, ny, nz) = size(z)
    x = [xi for xi in x, j in 1:ny, k in 1:nz]
    y = [yj for i in 1:nx, yj in y, k in 1:nz]

    # Get the time.
    t = round.(data["t"][:]; sigdigits) ./ time_unit_factor

    # Create the figure.
    figure = Figure()

    # Loop over outputs.
    row = 0
    for (variable, i, j, k, n) in fields
        row += 1
        column = 0

        # Round the time.
        tn = round(t[n]; sigdigits)

        # Get the label.
        label = LaTeXString(attrs(data[variable])["label"])

        if variable in ray_volume_properties
            # Get the ray-volume data.
            xr =
                round.(data["xr"][:, :, :, :, n]; sigdigits) ./
                space_unit_factor
            yr =
                round.(data["yr"][:, :, :, :, n]; sigdigits) ./
                space_unit_factor
            zr =
                round.(data["zr"][:, :, :, :, n]; sigdigits) ./
                space_unit_factor
            dxr =
                round.(data["dxr"][:, :, :, :, n]; sigdigits) ./
                space_unit_factor
            dyr =
                round.(data["dyr"][:, :, :, :, n]; sigdigits) ./
                space_unit_factor
            dzr =
                round.(data["dzr"][:, :, :, :, n]; sigdigits) ./
                space_unit_factor
            nr = round.(data["nr"][:, :, :, :, n]; sigdigits)
            phi = round.(data[variable][:, :, :, :, n]; sigdigits)

            # Plot in the x-y plane.
            if nx > 1 && ny > 1
                column += 2
                zk = round(sum(z[:, :, k]) / length(z[:, :, k]); sigdigits)
                add_scatter_plot!(
                    figure,
                    (;
                        colormap_name,
                        columns = (column - 1):column,
                        dx = dxr[:, :, :, k],
                        dy = dyr[:, :, :, k],
                        label,
                        mask = nr[:, :, :, k] .!= 0,
                        number,
                        phi = phi[:, :, :, k],
                        row,
                        title = L"t\approx%$tn\ \mathrm{%$time_unit},\quad z\approx%$zk\ \mathrm{%$space_unit}",
                        x = xr[:, :, :, k],
                        x_label = L"x_r\ [\mathrm{%$space_unit}]",
                        xmax = maximum(x),
                        xmin = minimum(x),
                        y = yr[:, :, :, k],
                        y_label = L"y_r\ [\mathrm{%$space_unit}]",
                        ymax = maximum(y),
                        ymin = minimum(y),
                    ),
                )
            end

            # Plot in the x-z plane.
            if nx > 1 && nz > 1
                column += 2
                yj = round(sum(y[:, j, :]) / length(y[:, j, :]); sigdigits)
                add_scatter_plot!(
                    figure,
                    (;
                        colormap_name,
                        columns = (column - 1):column,
                        dx = dxr[:, :, j, :],
                        dy = dzr[:, :, j, :],
                        label,
                        mask = nr[:, :, j, :] .!= 0,
                        number,
                        phi = phi[:, :, j, :],
                        row,
                        title = L"t\approx%$tn\ \mathrm{%$time_unit},\quad y\approx%$yj\ \mathrm{%$space_unit}",
                        x = xr[:, :, j, :],
                        x_label = L"x_r\ [\mathrm{%$space_unit}]",
                        xmax = maximum(x),
                        xmin = minimum(x),
                        y = zr[:, :, j, :],
                        y_label = L"z_r\ [\mathrm{%$space_unit}]",
                        ymax = maximum(z),
                        ymin = minimum(z),
                    ),
                )
            end

            # Plot in the y-z plane.
            if ny > 1 && nz > 1
                column += 2
                xi = round(sum(x[i, :, :]) / length(x[i, :, :]); sigdigits)
                add_scatter_plot!(
                    figure,
                    (;
                        colormap_name,
                        columns = (column - 1):column,
                        dx = dyr[:, i, :, :],
                        dy = dzr[:, i, :, :],
                        label,
                        mask = nr[:, i, :, :] .!= 0,
                        number,
                        phi = phi[:, i, :, :],
                        row,
                        title = L"t\approx%$tn\ \mathrm{%$time_unit},\quad x\approx%$xi\ \mathrm{%$space_unit}",
                        x = yr[:, i, :, :],
                        x_label = L"x_r\ [\mathrm{%$space_unit}]",
                        xmax = maximum(y),
                        xmin = minimum(y),
                        y = zr[:, i, :, :],
                        y_label = L"z_r\ [\mathrm{%$space_unit}]",
                        ymax = maximum(z),
                        ymin = minimum(z),
                    ),
                )
            end
        else
            # Get the variable.
            phi = round.(data[variable][:, :, :, n]; sigdigits)

            # Plot in the x-y plane.
            if nx > 1 && ny > 1
                column += 2
                zk = round(sum(z[:, :, k]) / length(z[:, :, k]); sigdigits)
                add_contour_plot!(
                    figure,
                    (;
                        background_color = :black,
                        colormap_name,
                        columns = (column - 1):column,
                        label,
                        number,
                        phi = phi[:, :, k],
                        row,
                        title = L"t\approx%$tn\ \mathrm{%$time_unit},\quad z\approx%$zk\ \mathrm{%$space_unit}",
                        x = x[:, :, k],
                        x_label = L"x\ [\mathrm{%$space_unit}]",
                        y = y[:, :, k],
                        y_label = L"y\ [\mathrm{%$space_unit}]",
                    ),
                )
            end

            # Plot in the x-z plane.
            if nx > 1 && nz > 1
                column += 2
                yj = round(sum(y[:, j, :]) / length(y[:, j, :]); sigdigits)
                add_contour_plot!(
                    figure,
                    (;
                        background_color = :black,
                        colormap_name,
                        columns = (column - 1):column,
                        label,
                        number,
                        phi = phi[:, j, :],
                        row,
                        title = L"t\approx%$tn\ \mathrm{%$time_unit},\quad y\approx%$yj\ \mathrm{%$space_unit}",
                        x = x[:, j, :],
                        x_label = L"x\ [\mathrm{%$space_unit}]",
                        y = z[:, j, :],
                        y_label = L"z\ [\mathrm{%$space_unit}]",
                    ),
                )
            end

            # Plot in the y-z plane.
            if ny > 1 && nz > 1
                column += 2
                xi = round(sum(x[i, :, :]) / length(x[i, :, :]); sigdigits)
                add_contour_plot!(
                    figure,
                    (;
                        background_color = :black,
                        colormap_name,
                        columns = (column - 1):column,
                        label,
                        number,
                        phi = phi[i, :, :],
                        row,
                        title = L"t\approx%$tn\ \mathrm{%$time_unit},\quad x\approx%$xi\ \mathrm{%$space_unit}",
                        x = y[i, :, :],
                        x_label = L"y\ [\mathrm{%$space_unit}]",
                        y = z[i, :, :],
                        y_label = L"z\ [\mathrm{%$space_unit}]",
                    ),
                )
            end
        end
    end

    # Resize, display and save the figure.
    resize_to_layout!(figure)
    display(figure)
    save(file, figure)

    return
end
