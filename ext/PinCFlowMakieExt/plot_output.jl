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
    animate::Bool = false,
    colormap_name::Symbol = :seismic,
    color_tick_format::AbstractString = "{:.1E}",
    display_figure::Bool = true,
    framerate::Real = 1,
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

  - `data_file`: HDF5 file with PinCFlow.jl output data.

  - `fields`: Either tuples of a variable name and a temporal index (for 2D data), or tuples of a variable name, three fractions which define the relative positions of the ``\\hat{y}``-``\\hat{z}``, ``\\hat{x}``-``\\hat{z}``, and ``\\hat{x}``-``\\hat{y}`` planes, and a temporal index (for 3D data).

# Keywords

  - `animate`: Switch for creating an animation instead of a multi-row figure.

  - `colormap_name`: Colormap of choice.

  - `color_tick_format`: Format string for the ticks of the colorbars.

  - `display_figure`: Switch for showing the figure with Makie.jl's `display` function.

  - `framerate`: Number of frames per second in the animation.

  - `number`: Number of contour levels.

  - `space_unit`: Unit used for the coordinates. Must be `:km` or `:m`.

  - `time_unit`: Unit used for the time. Must be `:d`, `:h`, `:min` or `:s`.

  - `x_tick_format`: Format string for the ticks of the ``x``-axis.

  - `y_tick_format`: Format string for the ticks of the ``y``-axis.

  - `z_tick_format`: Format string for the ticks of the ``z``-axis.

# See also

  - [`PinCFlow.set_visualization_theme!`](@ref)

  - [`PinCFlowMakieExt.plot_snapshot!`](@ref)
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
    animate::Bool = false,
    colormap_name::Symbol = :seismic,
    color_tick_format::AbstractString = "{:.1E}",
    display_figure::Bool = true,
    framerate::Real = 1,
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
    h5open(data_file) do data

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

        # Create the NamedTuple that is to be passed to plot_snapshot!.
        input = (;
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
        )

        # Create an animation or plot snapshots.
        if animate
            plot_snapshot!(figure, fields[1], input, 1)
            resize_to_layout!(figure)
            record(figure, plot_file, fields; framerate) do field
                plot_snapshot!(figure, field, input, 1)
                return
            end
        else
            row = 0
            for field in fields
                row += 1
                plot_snapshot!(figure, field, input, row)
            end
            resize_to_layout!(figure)
            display_figure && display(figure)
            save(plot_file, figure)
        end

        return
    end

    return
end
