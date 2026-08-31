# src/Examples/periodic_hill.jl

"""
```julia
periodic_hill(;
    display_figure::Bool = true,
    npx::Integer = 1,
    npy::Integer = 1,
    npz::Integer = 1,
    output_file::AbstractString = "periodic_hill.h5",
    plot_file::AbstractString = "periodic_hill.svg",
    prepare_restart::Bool = false,
    visualize::Bool = true,
    x_size::Integer = 20,
    y_size::Integer = 20,
    z_size::Integer = 20,
)
```

Run the periodic-hill example simulation.

# Keywords

The keywords are analogous to those of `cold_bubble`.
"""
function periodic_hill end

function periodic_hill(;
    display_figure::Bool = true,
    npx::Integer = 1,
    npy::Integer = 1,
    npz::Integer = 1,
    output_file::AbstractString = "periodic_hill.h5",
    plot_file::AbstractString = "periodic_hill.svg",
    prepare_restart::Bool = false,
    visualize::Bool = true,
    x_size::Integer = 20,
    y_size::Integer = 1,
    z_size::Integer = 20,
)
    h0 = 500
    l0 = 10000

    lz = 20000
    zr = 10000

    atmosphere = AtmosphereNamelist(;
        background = :StableStratification,
        coriolis_frequency = 0.0,
        initial_u = (x, y, z) -> 10.0,
        model = :Boussinesq,
    )

    domain = DomainNamelist(; lx = 20000, lz, npx, npz, x_size, y_size, z_size)

    grid = GridNamelist(;
        resolved_topography = (x, y) -> h0 / 2 * (1 + cos(pi / l0 * x)),
    )

    output =
        OutputNamelist(; output_file, output_variables = [:w], prepare_restart)

    sponge = SpongeNamelist(;
        rhs_sponge = (x, y, z, t, dt) ->
            z >= zr ? sin(pi / 2 * (z - zr) / (lz - zr))^2 / dt : 0.0,
    )

    integrate(Namelists(; atmosphere, domain, grid, output, sponge))

    if visualize && MPI.Comm_rank(MPI.COMM_WORLD) == 0
        plot_output(plot_file, output_file, (:w, 2); display_figure)
    end

    return
end
