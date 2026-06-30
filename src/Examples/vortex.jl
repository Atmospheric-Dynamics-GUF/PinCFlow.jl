# src/Examples/vortex.jl

function vortex(;
    display_figure::Bool = true,
    npx::Integer = 1,
    npy::Integer = 1,
    output_file::AbstractString = "vortex.h5",
    plot_file::AbstractString = "vortex.svg",
    prepare_restart::Bool = false,
    visualize::Bool = true,
    x_size::Integer = 20,
    y_size::Integer = 20,
    z_size::Integer = 1,
)
    lx = 20000
    ly = 20000

    rx = lx / 4
    ry = ly / 4

    atmosphere = AtmosphereNamelist(;
        model = :Boussinesq,
        background = :NeutralStratification,
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

    domain = DomainNamelist(; lx, ly, npx, npy, x_size, y_size, z_size)

    output = OutputNamelist(;
        output_file,
        output_variables = [:chi],
        prepare_restart,
    )

    tracer = TracerNamelist(;
        tracer_setup = :TracerOn,
        initial_chi = (x, y, z) -> begin
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
        plot_output(plot_file, output_file, (:chi, 2); display_figure)
    end

    return
end
