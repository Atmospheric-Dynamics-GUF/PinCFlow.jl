# src/Examples/hot_bubble.jl

"""
```julia
hot_bubble(;
    display_figure::Bool = true,
    npx::Integer = 1,
    npy::Integer = 1,
    npz::Integer = 1,
    output_file::AbstractString = "hot_bubble.h5",
    plot_file::AbstractString = "hot_bubble.svg",
    prepare_restart::Bool = false,
    visualize::Bool = true,
    x_size::Integer = 20,
    y_size::Integer = 20,
    z_size::Integer = 20,
)
```

Run the hot-bubble example simulation.

# Keywords

The keywords are analogous to those of `cold_bubble`.
"""
function hot_bubble end

function hot_bubble(;
    display_figure::Bool = true,
    npx::Integer = 1,
    npy::Integer = 1,
    npz::Integer = 1,
    output_file::AbstractString = "hot_bubble.h5",
    plot_file::AbstractString = "hot_bubble.svg",
    prepare_restart::Bool = false,
    visualize::Bool = true,
    x_size::Integer = 20,
    y_size::Integer = 20,
    z_size::Integer = 20,
)
    lx = 10000
    ly = 10000
    lz = 10000

    rx = lx / 4
    ry = ly / 4
    rz = lz / 4

    atmosphere = AtmosphereNamelist(;
        background = :Isentropic,
        initial_rhop = (x, y, z) -> begin
            r = sqrt((x / rx)^2 + (y / ry)^2 + ((z - 3 * rz) / rz)^2)
            if r <= 1
                return -0.005 * (1 + cos(pi * r))
            else
                return 0.0
            end
        end,
        model = :Compressible,
    )

    discretization = DiscretizationNamelist(; dtmax = 60)

    domain = DomainNamelist(; lx, ly, lz, npx, npy, npz, x_size, y_size, z_size)

    output = OutputNamelist(;
        output_file,
        output_interval = 300,
        output_variables = [:thetap],
        prepare_restart,
        tmax = 300,
    )

    integrate(Namelists(; atmosphere, discretization, domain, output))

    if visualize && MPI.Comm_rank(MPI.COMM_WORLD) == 0
        plot_output(
            plot_file,
            output_file,
            (:thetap, 0.5, 0.5, 0.75, 2);
            display_figure,
            time_unit = :min,
        )
    end

    return
end
