# PinCFlowExamples.jl/src/vortex.jl

function vortex(;
    x_size::Int64 = 40,
    y_size::Int64 = 40,
    npx::Int64 = 3,
    npy::Int64 = 3,
    output_file::String = "vortex.h5",
    prepare_restart::Bool = false,
    output_steps::Bool = false,
    visualize::Bool = true,
)
    lx = 20000.0
    ly = 20000.0

    rx = lx / 4
    ry = ly / 4

    atmosphere = AtmosphereNamelist(;
        model = Boussinesq(),
        background = NeutralStratification(),
        initial_u = (x, y, z) -> begin
            r = sqrt((x / rx)^2 + (y / ry)^2)
            if r <= 1
                return -5 * y / ry * (1 + cos(pi * r)) / 2
            else
                return 0.0
            end
        end,
        initial_v = (x, y, z) -> begin
            r = sqrt((x / rx)^2 + (y / ry)^2)
            if r <= 1
                return 5 * x / rx * (1 + cos(pi * r)) / 2
            else
                return 0.0
            end
        end,
    )

    domain = DomainNamelist(; x_size, y_size, lx, ly, npx, npy)

    output = OutputNamelist(;
        output_file,
        output_variables = (:chi,),
        prepare_restart,
        output_steps,
    )

    tracer = TracerNamelist(;
        tracer_setup = TracerOn(),
        initial_tracer = (x, y, z) -> begin
            r = sqrt(((abs(x) - rx) / rx)^2 + (y / ry)^2)
            if r <= 1
                return sign(x) * (1 + cos(pi * r)) / 2
            else
                return 0.0
            end
        end,
    )

    integrate(Namelists(; atmosphere, domain, output, tracer))

    if visualize && MPI.Comm_rank(MPI.COMM_WORLD) == 0
        h5open(output_file) do data
            plot_output(
                "PinCFlowExamples.jl/results/vortex.svg",
                data,
                ("chi", 1, 1, 1, 2);
            )
            return
        end
    end

    return
end
