"""
```julia
plot_snapshot!(
    figure::Figure,
    field::Union{
        Tuple{Symbol, <:Integer},
        Tuple{Symbol, <:Real, <:Real, <:Real, <:Integer},
    },
    input::NamedTuple,
    row::Integer,
)
```

Plot the snapshot specified by `field` in the row `row` of the figure `figure`.

# Arguments

  - `figure`: Figure to modify.

  - `field`: Field to plot.

  - `input`: Data and specifications for the plot.

  - `row`: Row of `figure` to add the plot to.

# See also

  - [`PinCFlowMakieExt.add_scatter_plot!`](@ref)

  - [`PinCFlowMakieExt.add_contour_plot!`](@ref)
"""
function plot_snapshot! end

@ivy function plot_snapshot!(
    figure::Figure,
    field::Union{
        Tuple{Symbol, <:Integer},
        Tuple{Symbol, <:Real, <:Real, <:Real, <:Integer},
    },
    input::NamedTuple,
    row::Integer,
)
    (;
        colormap_name,
        color_tick_format,
        data,
        digits,
        number,
        nx,
        ny,
        nz,
        plane_in_title,
        ray_volume_properties,
        space_unit,
        space_unit_factor,
        t,
        time_unit,
        x,
        x_tick_format,
        y,
        y_tick_format,
        z,
        z_tick_format,
    ) = input

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
    label = LaTeXString(attrs(data[string(variable)])["label"])

    column = 0

    if variable in ray_volume_properties
        # Get the ray-volume data.
        xr = data["xr"][:, :, :, :, n] ./ space_unit_factor
        yr = data["yr"][:, :, :, :, n] ./ space_unit_factor
        zr = data["zr"][:, :, :, :, n] ./ space_unit_factor
        dxr = data["dxr"][:, :, :, :, n] ./ space_unit_factor
        dyr = data["dyr"][:, :, :, :, n] ./ space_unit_factor
        dzr = data["dzr"][:, :, :, :, n] ./ space_unit_factor
        nr = data["nr"][:, :, :, :, n]
        phi = data[string(variable)][:, :, :, :, n]

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
        phi = data[string(variable)][:, :, :, n]

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

    return
end
