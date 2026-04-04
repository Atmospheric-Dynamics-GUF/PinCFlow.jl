# PinCFlowExamples.jl/src/wave_packet.jl

function wave_packet(;
    x_size::Integer = 40,
    y_size::Integer = 40,
    z_size::Integer = 80,
    npx::Integer = 3,
    npy::Integer = 3,
    npz::Integer = 3,
    output_file::AbstractString = "wave_packet.h5",
    prepare_restart::Bool = false,
    output_steps::Bool = false,
    visualize::Bool = true,
)
    lx = 20000.0
    ly = 20000.0
    lz = 40000.0

    parameters = (
        k = 16 * pi / lx,
        l = 16 * pi / ly,
        m = 32 * pi / lz,
        rx = 0.25,
        ry = 0.25,
        rz = 0.25,
        x0 = 0.0,
        y0 = 0.0,
        z0 = 20000.0,
        a0 = 0.05,
    )

    background = Realistic()
    coriolis_frequency = 0.0001

    atmosphere = AtmosphereNamelist(; background, coriolis_frequency)

    domain = DomainNamelist(; x_size, y_size, z_size, lx, ly, lz, npx, npy, npz)

    state = State(Namelists(; atmosphere, domain))
    (; g) = state.constants

    atmosphere = AtmosphereNamelist(;
        background,
        coriolis_frequency,
        initial_u = (x, y, z) -> real(
            uhat(state, parameters, x, y, z) *
            exp(1im * phi(parameters, x, y, z)),
        ),
        initial_v = (x, y, z) -> real(
            vhat(state, parameters, x, y, z) *
            exp(1im * phi(parameters, x, y, z)),
        ),
        initial_w = (x, y, z) -> real(
            what(state, parameters, x, y, z) *
            exp(1im * phi(parameters, x, y, z)),
        ),
        initial_pip = (x, y, z) -> real(
            pihat(state, parameters, x, y, z) *
            exp(1im * phi(parameters, x, y, z)),
        ),
        initial_thetap = (x, y, z) ->
            real(
                bhat(state, parameters, x, y, z) *
                exp(1im * phi(parameters, x, y, z)),
            ) / g * thetabar(state, x, y, z),
        buoyancy_initialization = :initial_thetap,
    )

    output = OutputNamelist(;
        output_file,
        output_variables = (:u, :v, :w),
        prepare_restart,
        output_steps,
        output_interval = 900,
        tmax = 900,
    )

    integrate(Namelists(; atmosphere, domain, output))

    if visualize && MPI.Comm_rank(MPI.COMM_WORLD) == 0
        h5open(output_file) do data
            plot_output(
                "PinCFlowExamples.jl/results/wave_packet.svg",
                data,
                ("u", 20, 20, 40, 2),
                ("v", 20, 20, 40, 2),
                ("w", 20, 20, 40, 2);
                time_unit = "min",
            )
            return
        end
    end

    return
end
