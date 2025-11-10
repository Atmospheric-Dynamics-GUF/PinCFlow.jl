function plot_contours(
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

    # Define the labels.
    labels = Dict(
        "rhobar" => L"\overline{\rho}\,[\mathrm{kg\,m^{-3}}]",
        "thetabar" => L"\overline{\theta}\,[\mathrm{K}]",
        "n2" => L"N^2\,[\mathrm{s^{-2}}]",
        "p" => L"P\,[\mathrm{kg\,m^{-3}\,K}]",
        "rhop" => L"\rho'\,[\mathrm{kg\,m^{-3}}]",
        "u" => L"u\,[\mathrm{m\,s^{-1}}]",
        "us" => L"u_\mathrm{s}\,[\mathrm{m\,s^{-1}}]",
        "v" => L"v\,[\mathrm{m\,s^{-1}}]",
        "vs" => L"v_\mathrm{s}\,[\mathrm{m\,s^{-1}}]",
        "w" => L"w\,[\mathrm{m\,s^{-1}}]",
        "ws" => L"w_\mathrm{s}\,[\mathrm{m\,s^{-1}}]",
        "wt" => L"\widehat{w}\,[\mathrm{m\,s^{-1}}]",
        "wts" => L"\widehat{w}_\mathrm{s}\,[\mathrm{m\,s^{-1}}]",
        "thetap" => L"\theta'\,[\mathrm{K}]",
        "pip" => L"\pi'",
        "chi" => L"\chi",
        "dudt" =>
            L"\left(\frac{\partial \rho u_\mathrm{b}}{\partial t}\right)_\mathrm{w}\,[\mathrm{kg\,m^{-2}\,s^{-2}}]",
        "dvdt" =>
            L"\left(\frac{\partial \rho v_\mathrm{b}}{\partial t}\right)_\mathrm{w}\,[\mathrm{kg\,m^{-2}\,s^{-2}}]",
        "dthetadt" =>
            L"\left(\frac{\partial P_\mathrm{b}}{\partial t}\right)_\mathrm{w}\,[\mathrm{kg\,m^{-3}\,K\,s^{-1}}]",
        "dchidt" =>
            L"\left(\frac{\partial \rho \chi_\mathrm{b}}{\partial t}\right)_\mathrm{w}\,[\mathrm{kg\,m^{-3}\,s^{-1}}]",
        "uchi" =>
            L"\overline{\rho}\langle u'\chi' \rangle\,[\mathrm{kg\,m^{-2}\,s^{-1}}]",
        "vchi" =>
            L"\overline{\rho}\langle v'\chi' \rangle\,[\mathrm{kg\,m^{-2}\,s^{-1}}]",
        "wchi" =>
            L"\overline{\rho}\langle w'\chi' \rangle\,[\mathrm{kg\,m^{-2}\,s^{-1}}]",
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

        # Get the label.
        label = labels[variable]

        # Get the variable.
        phi = data[variable][:, :, :, :]

        # Round the time.
        tn = round(t[n]; digits = 1)

        # Plot in the x-y plane.
        if nx > 1 && ny > 1
            column += 2
            @ivy zk = round(sum(z[:, :, k]) / length(z[:, :, k]); digits = 1)
            axis = Axis(
                figure[row, column - 1];
                title = L"t\approx%$tn\,\mathrm{%$time_unit},\quad z\approx%$zk\,\mathrm{%$space_unit}",
                xlabel = L"x\,[\mathrm{%$space_unit}]",
                ylabel = L"y\,[\mathrm{%$space_unit}]",
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
            @ivy yj = round(sum(y[:, j, :]) / length(y[:, j, :]); digits = 1)
            axis = Axis(
                figure[row, column - 1];
                title = L"t\approx%$tn\,\mathrm{%$time_unit},\quad y\approx%$yj\,\mathrm{%$space_unit}",
                xlabel = L"x\,[\mathrm{%$space_unit}]",
                ylabel = L"z\,[\mathrm{%$space_unit}]",
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
            @ivy xi = round(sum(x[i, :, :]) / length(x[i, :, :]); digits = 1)
            axis = Axis(
                figure[row, column - 1];
                title = L"t\approx%$tn\,\mathrm{%$time_unit},\quad x\approx%$xi\,\mathrm{%$space_unit}",
                xlabel = L"y\,[\mathrm{%$space_unit}]",
                ylabel = L"z\,[\mathrm{%$space_unit}]",
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
                ticks = levels,
                tickformat = "{:9.2E}",
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
