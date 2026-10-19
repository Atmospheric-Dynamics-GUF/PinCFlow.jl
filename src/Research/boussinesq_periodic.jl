# src/Examples/boussinesq_periodic.jl

function boussinesq_periodic(;
    display_figure::Bool = true,
    npx::Integer = 1,
    npy::Integer = 1,
    npz::Integer = 1,
    output_file::AbstractString = "boussinesq_periodic.h5",
    plot_file::AbstractString = "boussinesq_periodic.svg",
    prepare_restart::Bool = false,
    visualize::Bool = true,
    x_size::Integer = 216,
    y_size::Integer = 16,
    z_size::Integer = 320,
)
    lx = 400.0e3
    ly = 30.0e3
    lz = 20.0e3

    parameters = (
        k = 0.0,
        l = 2 * pi / 30.0e3,
        m = 2 * pi / 1.0e3,
        rx = 0.0,
        ry = 0.0,
        rz = 0.0,
        x0 = 0.0,
        y0 = 0.0,
        z0 = 0.0,
        a0 = 0.7,
        lz = lz,
    )
    (; m) = parameters

    model = :Boussinesq
    background = :StableStratification
    coriolis_frequency = 2.0e-3
    buoyancy_frequency = 2.0e-2

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
        vertical_boundary_condition = :Periodic,
    )

    atmosphere = AtmosphereNamelist(;
        model,
        background,
        coriolis_frequency,
        buoyancy_frequency,
    )

    state = State(Namelists(; atmosphere, domain))
    (; g) = state.constants

    atmosphere = AtmosphereNamelist(;
        model,
        background,
        coriolis_frequency,
        buoyancy_frequency,
        buoyancy_initialization = :initial_thetap,
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
        output_interval = 1.0,
        output_variables = [:u, :v, :w, :tke, :rhop, :tracerdiffusion],
        prepare_restart,
        tmax = 1.0,
    )

    function chi_largescale(x::Real, z::Real)
        return 0.25 *
               (1 + cos(2 * pi * (z - lz / 2) / lz)) *
               (1 + cos(2 * pi * x / lx))
    end

    function dchidx(x::Real, z::Real)
        return -pi / (2 * lx) *
               (1 + cos(2 * pi * (z - lz / 2) / lz)) *
               sin(2 * pi * x / lx)
    end

    function dchidz(x::Real, z::Real)
        return -pi / (2 * lz) *
               sin(2 * pi * (z - lz / 2) / lz) *
               (1 + cos(2 * pi * x / lx))
    end

    function chihat(
        state::State,
        parameters::NamedTuple,
        x::Real,
        y::Real,
        z::Real,
    )
        return -1im / omega(state, parameters, x, y, z) * (
            uhat(state, parameters, x, y, z) * dchidx(x, z) +
            what(state, parameters, x, y, z) * dchidz(x, z)
        )
    end

    tracer = TracerNamelist(;
        tracer_setup = :TracerOn,
        initial_chi = (x, y, z) ->
            chi_largescale(x, z) + real(
                chihat(state, parameters, x, y, z) *
                exp(1im * phi(parameters, x, y, z)),
            ),
    )

    lt = 30.0

    turbulence = TurbulenceNamelist(;
        initial_tke = (x, y, z) ->
            max(
                0,
                lt * (
                    lt * m^2 / 2 * (
                        abs(uhat(state, parameters, x, y, z))^2 +
                        abs(vhat(state, parameters, x, y, z))^2 - real(
                            (
                                uhat(state, parameters, x, y, z) *
                                uhat(state, parameters, x, y, z) +
                                vhat(state, parameters, x, y, z) *
                                vhat(state, parameters, x, y, z)
                            ) * exp(2im * phi(parameters, x, y, z)),
                        )
                    ) -
                    lt * (
                        n2(state, x, y, z) + real(
                            1im *
                            m *
                            bhat(state, parameters, x, y, z) *
                            exp(1im * phi(parameters, x, y, z)),
                        )
                    )
                ),
            ) / 2,
    )

    integrate(Namelists(; atmosphere, domain, output, tracer, turbulence))

    if visualize && MPI.Comm_rank(MPI.COMM_WORLD) == 0
        plot_output(
            plot_file,
            output_file,
            (:u, 0.5, 0.5, 0.5, 2),
            (:v, 0.5, 0.5, 0.5, 2),
            (:w, 0.5, 0.5, 0.5, 2);
            display_figure,
            time_unit = :min,
        )
    end

    return
end
