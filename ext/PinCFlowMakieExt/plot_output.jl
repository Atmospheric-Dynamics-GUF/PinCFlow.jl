@ivy function plot_output(
    plot_file::AbstractString,
    data_file::AbstractString,
    fields::Vararg{
        Union{
            Tuple{Symbol, <:Integer},
            Tuple{Symbol, <:Real, <:Real, <:Real, <:Integer},
        },
    };
    colormap_name::Symbol = :seismic,
    color_tick_format::AbstractString = "{:.2E}",
    number::Integer = 10,
    significant_digits::Integer = 3,
    space_unit::Symbol = :km,
    time_unit::Symbol = :h,
    x_tick_format::AbstractString = "{:.0f}",
    y_tick_format::AbstractString = "{:.0f}",
    z_tick_format::AbstractString = "{:.0f}",
)
    set_visualization_theme!()

    # Store the ray-volume property names.
    ray_volume_properties =
        (:xr, :yr, :zr, :dxr, :dyr, :dzr, :kr, :lr, :mr, :dkr, :dlr, :dmr, :nr)

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

    # Create the figure.
    figure = h5open(data_file) do data

        # Set the grid.
        x = round.(data["x"][:]; sigdigits) ./ space_unit_factor
        y = round.(data["y"][:]; sigdigits) ./ space_unit_factor
        z = round.(data["z"][:, :, :]; sigdigits) ./ space_unit_factor
        (nx, ny, nz) = size(z)
        x = [xi for xi in x, j in 1:ny, k in 1:nz]
        y = [yj for i in 1:nx, yj in y, k in 1:nz]

        # Get the time.
        t = round.(data["t"][:]; sigdigits) ./ time_unit_factor

        figure = Figure()

        # Loop over outputs.
        row = 0
        for field in fields
            row += 1
            column = 0

            # Check if the fields specification is correct.
            if nx > 1 && ny > 1 && nz > 1 && length(field) != 5
                error("Incorrect fields specification for 3D data!")
            end
            if length(field) == 5
                for index in 2:4
                    if field[index] < 0 || field[index] > 1
                        error("Incorrect plane specification!")
                    end
                end
            end

            # Determine the data slices.
            (variable, n) = (field[1], field[end])
            if length(field) == 5
                (i, j, k) = ceil.(Int64, field[2:4] .* (nx, ny, nz))
            else
                i = j = k = 1
            end

            # Round the time.
            tn = round(t[n]; sigdigits)

            # Get the label.
            label = LaTeXString(attrs(data[string(variable)])["label"])

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
                phi =
                    round.(data[string(variable)][:, :, :, :, n]; sigdigits)

                # Plot in the x-y plane.
                if nx > 1 && ny > 1
                    column += 2
                    zk = round(sum(z[:, :, k]) / length(z[:, :, k]); sigdigits)
                    add_scatter_plot!(
                        figure,
                        (;
                            colormap_name,
                            color_tick_format,
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
                            x_tick_format,
                            y = yr[:, :, :, k],
                            y_label = L"y_r\ [\mathrm{%$space_unit}]",
                            ymax = maximum(y),
                            ymin = minimum(y),
                            y_tick_format,
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
                            color_tick_format,
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
                            x_tick_format,
                            y = zr[:, :, j, :],
                            y_label = L"z_r\ [\mathrm{%$space_unit}]",
                            ymax = maximum(z),
                            ymin = minimum(z),
                            y_tick_format = z_tick_format,
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
                            color_tick_format,
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
                            x_tick_format = y_tick_format,
                            y = zr[:, i, :, :],
                            y_label = L"z_r\ [\mathrm{%$space_unit}]",
                            ymax = maximum(z),
                            ymin = minimum(z),
                            y_tick_format = z_tick_format,
                        ),
                    )
                end
            else
                # Get the variable.
                phi = round.(data[string(variable)][:, :, :, n]; sigdigits)

                # Plot in the x-y plane.
                if nx > 1 && ny > 1
                    column += 2
                    zk = round(sum(z[:, :, k]) / length(z[:, :, k]); sigdigits)
                    add_contour_plot!(
                        figure,
                        (;
                            background_color = :black,
                            colormap_name,
                            color_tick_format,
                            columns = (column - 1):column,
                            label,
                            number,
                            phi = phi[:, :, k],
                            row,
                            title = L"t\approx%$tn\ \mathrm{%$time_unit},\quad z\approx%$zk\ \mathrm{%$space_unit}",
                            x = x[:, :, k],
                            x_label = L"x\ [\mathrm{%$space_unit}]",
                            x_tick_format,
                            y = y[:, :, k],
                            y_label = L"y\ [\mathrm{%$space_unit}]",
                            y_tick_format,
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
                            color_tick_format,
                            columns = (column - 1):column,
                            label,
                            number,
                            phi = phi[:, j, :],
                            row,
                            title = L"t\approx%$tn\ \mathrm{%$time_unit},\quad y\approx%$yj\ \mathrm{%$space_unit}",
                            x = x[:, j, :],
                            x_label = L"x\ [\mathrm{%$space_unit}]",
                            x_tick_format,
                            y = z[:, j, :],
                            y_label = L"z\ [\mathrm{%$space_unit}]",
                            y_tick_format = z_tick_format,
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
                            color_tick_format,
                            columns = (column - 1):column,
                            label,
                            number,
                            phi = phi[i, :, :],
                            row,
                            title = L"t\approx%$tn\ \mathrm{%$time_unit},\quad x\approx%$xi\ \mathrm{%$space_unit}",
                            x = y[i, :, :],
                            x_label = L"y\ [\mathrm{%$space_unit}]",
                            x_tick_format = y_tick_format,
                            y = z[i, :, :],
                            y_label = L"z\ [\mathrm{%$space_unit}]",
                            y_tick_format = z_tick_format,
                        ),
                    )
                end
            end
        end

        return figure
    end

    # Resize, display and save the figure.
    resize_to_layout!(figure)
    display(figure)
    save(plot_file, figure)

    return
end
