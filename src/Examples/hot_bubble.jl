# src/Examples/hot_bubble.jl

function hot_bubble(;
    x_size::Integer = 40,
    z_size::Integer = 40,
    npx::Integer = 1,
    npz::Integer = 1,
    output_file::AbstractString = "hot_bubble.h5",
    prepare_restart::Bool = false,
    output_steps::Bool = false,
    visualize::Bool = true,
    plot_file::AbstractString = "examples/results/hot_bubble.svg",
)
    lx = 20000.0
    lz = 20000.0

    rx = lx / 8
    rz = lz / 8

    atmosphere = AtmosphereNamelist(;
        model = :Compressible,
        background = :Isentropic,
        initial_rhop = (x, y, z) -> begin
            r = sqrt((x / rx)^2 + ((z - 5 * rz) / rz)^2)
            if r <= 1
                return -0.005 * (1 + cos(pi * r))
            else
                return 0.0
            end
        end,
    )

    discretization = DiscretizationNamelist(; dtmax = 60.0)

    domain = DomainNamelist(; x_size, z_size, lx, lz, npx, npz)

    output = OutputNamelist(;
        output_file,
        output_variables = [:thetap],
        prepare_restart,
        output_steps,
    )

    integrate(Namelists(; atmosphere, discretization, domain, output))

    if visualize && MPI.Comm_rank(MPI.COMM_WORLD) == 0
        h5open(output_file) do data
            plot_output(plot_file, data, ("thetap", 1, 1, 1, 2))
            return
        end
    end

    return
end
