"""
```julia
plot_output(
    plot_file::AbstractString,
    data_file::AbstractString,
    fields::Vararg{
        Union{
            Tuple{Symbol, <:Integer},
            Tuple{Symbol, <:Real, <:Real, <:Real, <:Integer},
        },
    };
    colormap_name::Symbol = :seismic,
    color_tick_format::AbstractString = "{:.1E}",
    number::Integer = 10,
    space_unit::Symbol = :km,
    time_unit::Symbol = :h,
    x_tick_format::AbstractString = "{:.1f}",
    y_tick_format::AbstractString = "{:.1f}",
    z_tick_format::AbstractString = "{:.1f}",
)
```

Create contour plots of the dataset `variable` in `data`, display it and save it to `file`.

# Arguments

  - `plot_file`: File to save the plots to.

  - `data_file`: HDF5 or NetCDF4 file with PinCFlow.jl output data.

  - `fields`: Either tuples of a variable name and a temporal index (for 2D data), or tuples of a variable name, three fractions which define the relative positions of the ``\\hat{y}``-``\\hat{z}``, ``\\hat{x}``-``\\hat{z}``, and ``\\hat{x}``-``\\hat{y}`` planes, and a temporal index (for 3D data).

# Keywords

  - `colormap_name`: Colormap of choice.

  - `color_tick_format`: Format string for the ticks of the colorbars.

  - `display_figure`: Switch for showing the figure with Makie.jl's `display` function.

  - `number`: Number of contour levels.

  - `space_unit`: Unit used for the coordinates. Must be `:km` or `:m`.

  - `time_unit`: Unit used for the time. Must be `:d`, `:h`, `:min` or `:s`.

  - `x_tick_format`: Format string for the ticks of the ``x``-axis.

  - `y_tick_format`: Format string for the ticks of the ``y``-axis.

  - `z_tick_format`: Format string for the ticks of the ``z``-axis.

# See also

  - [`PinCFlow.set_visualization_theme!`](@ref)

  - [`PinCFlowMakieExt.add_scatter_plot!`](@ref)

  - [`PinCFlowMakieExt.add_contour_plot!`](@ref)
"""
plot_output

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
    color_tick_format::AbstractString = "{:.1E}",
    display_figure::Bool = true,
    number::Integer = 10,
    space_unit::Symbol = :km,
    time_unit::Symbol = :h,
    x_tick_format::AbstractString = "{:.1f}",
    y_tick_format::AbstractString = "{:.1f}",
    z_tick_format::AbstractString = "{:.1f}",
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
        error("Unknown space unit!")
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
        error("Unknown time unit!")
    end

    # Set the digits for rounding the time and plane positions.
    digits = 1

    # Create the figure.
    figure = NCDataset(data_file, "r") do data

        # Set the grid.
        x = data["x"][:] ./ space_unit_factor
        y = data["y"][:] ./ space_unit_factor
        z = data["z"][:, :, :] ./ space_unit_factor
        (nx, ny, nz) = size(z)
        x = [xi for xi in x, j in 1:ny, k in 1:nz]
        y = [yj for i in 1:nx, yj in y, k in 1:nz]

        # Determine whether the plane should be included in the title.
        plane_in_title = (nx > 1 && ny > 1 && nz > 1)

        # Get the time.
        t = data["t"][:] ./ time_unit_factor

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
            tn = round(t[n]; digits)

            # Get the label.
            label = LaTeXString(data[string(variable)].attrib["label"][:])

            if variable in ray_volume_properties
                # Get the ray-volume data.
                xr = data["xr"][:, :, :, :, n] ./ space_unit_factor
                yr = data["yr"][:, :, :, :, n] ./ space_unit_factor
                zr = data["zr"][:, :, :, :, n] ./ space_unit_factor
                dxr = data["dxr"][:, :, :, :, n] ./ space_unit_factor
                dyr = data["dyr"][:, :, :, :, n] ./ space_unit_factor
                dzr = data["dzr"][:, :, :, :, n] ./ space_unit_factor
                nr = data["nr"][:, :, :, :, n]
                phi = Array(data[string(variable)][:, :, :, :, n])

                # Plot in the x-y plane.
                if nx > 1 && ny > 1
                    column += 2
                    if plane_in_title
                        zk = round(sum(z[:, :, k]) / length(z[:, :, k]); digits)
                        title =
                            L"t\approx%$tn\ \mathrm{%$time_unit},\quad z\approx%$zk\ \mathrm{%$space_unit}"
                    else
                        title = L"t\approx%$tn\ \mathrm{%$time_unit}"
                    end
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
                            title,
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
                    if plane_in_title
                        yj = round(sum(y[:, j, :]) / length(y[:, j, :]); digits)
                        title =
                            L"t\approx%$tn\ \mathrm{%$time_unit},\quad y\approx%$yj\ \mathrm{%$space_unit}"
                    else
                        title = L"t\approx%$tn\ \mathrm{%$time_unit}"
                    end
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
                            title,
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
                    if plane_in_title
                        xi = round(sum(x[i, :, :]) / length(x[i, :, :]); digits)
                        title =
                            L"t\approx%$tn\ \mathrm{%$time_unit},\quad x\approx%$xi\ \mathrm{%$space_unit}"
                    else
                        title = L"t\approx%$tn\ \mathrm{%$time_unit}"
                    end
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
                            title,
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
                phi = Array(data[string(variable)][:, :, :, n])

                # Plot in the x-y plane.
                if nx > 1 && ny > 1
                    column += 2
                    if plane_in_title
                        zk = round(sum(z[:, :, k]) / length(z[:, :, k]); digits)
                        title =
                            L"t\approx%$tn\ \mathrm{%$time_unit},\quad z\approx%$zk\ \mathrm{%$space_unit}"
                    else
                        title = L"t\approx%$tn\ \mathrm{%$time_unit}"
                    end
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
                            title,
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
                    if plane_in_title
                        yj = round(sum(y[:, j, :]) / length(y[:, j, :]); digits)
                        title =
                            L"t\approx%$tn\ \mathrm{%$time_unit},\quad y\approx%$yj\ \mathrm{%$space_unit}"
                    else
                        title = L"t\approx%$tn\ \mathrm{%$time_unit}"
                    end
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
                            title,
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
                    if plane_in_title
                        xi = round(sum(x[i, :, :]) / length(x[i, :, :]); digits)
                        title =
                            L"t\approx%$tn\ \mathrm{%$time_unit},\quad x\approx%$xi\ \mathrm{%$space_unit}"
                    else
                        title = L"t\approx%$tn\ \mathrm{%$time_unit}"
                    end
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
                            title,
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
    display_figure && display(figure)
    save(plot_file, figure)

    return
end
