# src/Examples/periodic_hill.jl

function periodic_hill(;
    x_size::Integer = 40,
    z_size::Integer = 40,
    npx::Integer = 1,
    npz::Integer = 1,
    output_file::String = "periodic_hill.h5",
    prepare_restart::Bool = false,
    output_steps::Bool = false,
    visualize::Bool = true,
)
    h0 = 500.0
    l0 = 10000.0

    lz = 20000.0
    zr = 10000.0

    atmosphere = AtmosphereNamelist(;
        model = :Boussinesq,
        background = :StableStratification,
        coriolis_frequency = 0.0,
        initial_u = (x, y, z) -> 10.0,
    )

    domain = DomainNamelist(; x_size, z_size, lx = 20000.0, lz, npx, npz)

    grid = GridNamelist(;
        resolved_topography = (x, y) -> h0 / 2 * (1 + cos(pi / l0 * x)),
    )

    output = OutputNamelist(;
        output_file,
        output_variables = [:w],
        prepare_restart,
        output_steps,
    )

    sponge = SpongeNamelist(;
        rhs_sponge = (x, y, z, t, dt) ->
            z >= zr ? sin(pi / 2 * (z - zr) / (lz - zr))^2 / dt : 0.0,
    )

    integrate(Namelists(; atmosphere, domain, grid, output, sponge))

    if visualize && MPI.Comm_rank(MPI.COMM_WORLD) == 0
        h5open(output_file) do data
            plot_output(
                "examples/results/periodic_hill.svg",
                data,
                ("w", 1, 1, 1, 2);
            )
            return
        end
    end

    return
end
