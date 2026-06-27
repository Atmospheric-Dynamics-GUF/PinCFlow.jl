# src/Examples/wave_packet.jl

function wave_packet(;
    npx::Integer = 1,
    npy::Integer = 1,
    npz::Integer = 1,
    output_file::AbstractString = "wave_packet.h5",
    plot_file::AbstractString = "wave_packet.svg",
    prepare_restart::Bool = false,
    visualize::Bool = true,
    x_size::Integer = 20,
    y_size::Integer = 20,
    z_size::Integer = 20,
)
    z0 = 10000

    lx = 20000
    ly = 20000
    lz = 30000

    parameters = (
        k = 8 * pi / lx,
        l = 8 * pi / ly,
        m = 8 * pi / (lz - z0),
        rx = 0.5,
        ry = 0.5,
        rz = 0.5,
        x0 = 0.0,
        y0 = 0.0,
        z0 = 20000,
        a0 = 0.05,
    )

    background = :Realistic
    coriolis_frequency = 0.0001

    atmosphere = AtmosphereNamelist(; background, coriolis_frequency)

    domain = DomainNamelist(;
        lx,
        ly,
        lz,
        npx,
        npy,
        npz,
        x_size,
        y_size,
        z_size,
    )

    grid = GridNamelist(; resolved_topography = (x, y) -> z0)

    state = State(Namelists(; atmosphere, domain, grid))
    (; g) = state.constants

    atmosphere = AtmosphereNamelist(;
        background,
        buoyancy_initialization = :initial_thetap,
        coriolis_frequency,
        initial_pip = (x, y, z) -> real(
            pihat(state, parameters, x, y, z) *
            exp(1im * phi(parameters, x, y, z)),
        ),
        initial_thetap = (x, y, z) ->
            real(
                bhat(state, parameters, x, y, z) *
                exp(1im * phi(parameters, x, y, z)),
            ) / g * thetabar(state, x, y, z),
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
    )

    output = OutputNamelist(;
        output_file,
        output_interval = 600,
        output_variables = [:u, :v, :w],
        prepare_restart,
        tmax = 600,
    )

    integrate(Namelists(; atmosphere, domain, grid, output))

    if visualize && MPI.Comm_rank(MPI.COMM_WORLD) == 0
        plot_output(
            plot_file,
            output_file,
            (:u, 0.5, 0.5, 0.5, 2),
            (:v, 0.5, 0.5, 0.5, 2),
            (:w, 0.5, 0.5, 0.5, 2);
            time_unit = :min,
        )
    end

    return
end
