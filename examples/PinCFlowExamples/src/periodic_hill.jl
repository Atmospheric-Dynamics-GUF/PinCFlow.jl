# examples/PinCFlowExamples/src/periodic_hill.jl

function periodic_hill(;
    x_size::Integer = 40,
    z_size::Integer = 40,
    npx::Integer = 3,
    npz::Integer = 3,
    output_file::AbstractString = "periodic_hill.h5",
    output_steps::Bool = false,
    output_variables::Tuple{Vararg{Symbol}} = (:w,),
    prepare_restart::Bool = false,
    visualize::Bool = true,
)
    h0 = 500.0
    l0 = 10000.0

    lz = 20000.0
    zr = 10000.0

    atmosphere = AtmosphereNamelist(;
        model = Boussinesq(),
        background = StableStratification(),
        coriolis_frequency = 0.0,
        initial_u = (x, y, z) -> 10.0,
    )

    domain = DomainNamelist(; x_size, z_size, lx = 20000.0, lz, npx, npz)

    grid = GridNamelist(;
        resolved_topography = (x, y) -> h0 / 2 * (1 + cos(pi / l0 * x)),
    )

    output = OutputNamelist(;
        output_file,
        output_steps,
        output_variables,
        prepare_restart,
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
