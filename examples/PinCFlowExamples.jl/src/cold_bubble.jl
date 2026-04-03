# examples/PinCFlowExamples.jl/src/cold_bubble.jl

function cold_bubble(;
    x_size::Integer = 40,
    z_size::Integer = 40,
    npx::Integer = 3,
    npz::Integer = 3,
    output::OutputNamelist = OutputNamelist(;
        output_file = "cold_bubble.h5",
        output_variables = (:thetap,),
    ),
    visualize::Bool = true,
)
    lx = 20000.0
    lz = 20000.0

    rx = lx / 8
    rz = lz / 8

    atmosphere = AtmosphereNamelist(;
        background = Isentropic(),
        initial_rhop = (x, y, z) -> begin
            r = sqrt((x / rx)^2 + ((z - 3 * rz) / rz)^2)
            if r <= 1
                return 0.005 * (1 + cos(pi * r))
            else
                return 0.0
            end
        end,
    )

    discretization = DiscretizationNamelist(; dtmax = 60.0)

    domain = DomainNamelist(; x_size, z_size, lx, lz, npx, npz)

    integrate(Namelists(; atmosphere, discretization, domain, output))

    if visualize && MPI.Comm_rank(MPI.COMM_WORLD) == 0
        h5open(output_file) do data
            plot_output(
                "examples/results/cold_bubble.svg",
                data,
                ("thetap", 1, 1, 1, 2);
            )
            return
        end
    end

    return
end
