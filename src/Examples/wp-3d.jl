# src/Examples/wp-3d.jl

function wp_3d(;
    x_size::Integer = 512,
    y_size::Integer = 16,
    z_size::Integer = 1060,
    npx::Integer = 1,
    npy::Integer = 1,
    npz::Integer = 1,
    output_file::AbstractString = "wp-3d.h5",
    visualize::Bool = false,
    prepare_restart::Bool = false,
)
    lx = 9000e3
    ly = 300e3
    lz = 100e3

    parameters = (
        k = 0,
        l = 2 * pi / 300e3,
        m = 2 * pi / 3e3,
        rx = 1500e3,
        ry = 0,
        rz = 5e3,
        x0 = 0.0,
        y0 = 0.0,
        z0 = 30e3,
        a0 = 1.0,
        threedim = true,
    )

    background = :Isothermal
    model = :Compressible
    coriolis_frequency = 1e-4

    atmosphere = AtmosphereNamelist(; background, coriolis_frequency, model)

    domain = DomainNamelist(; x_size, y_size, z_size, lx, ly, lz, npx, npy, npz)

    state = State(Namelists(; atmosphere, domain))
    (; g) = state.constants
    (; lturb) = state.turbulence.turbulenceconstants

    atmosphere = AtmosphereNamelist(;
        background,
        coriolis_frequency,
        model,
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

    domain = DomainNamelist(; x_size, y_size, z_size, lx, ly, lz, npx, npy, npz)

    output = OutputNamelist(;
        output_variables = [:u, :v, :w, :rhop],
        output_file,
        prepare_restart,
        tmax = 100.0,
        output_interval = 100.0,
    )

    discretization = DiscretizationNamelist(; dtmax = 100)

    turbulence = TurbulenceNamelist(;
        turbulence_scheme = :NoTurbulence,
        initial_tke = (x, y, z) -> max(
            10.e-5,
            real(
                lturb^2.0 * (
                    m^2 / 2 * (
                        abs(uhat(state, parameters, x, y, z))^2 +
                        abs(vhat(state, parameters, x, y, z))^2 - real(
                            (
                                uhat(state, parameters, x, y, z)^2 +
                                vhat(state, parameters, x, y, z)^2
                            ) *
                            exp(2im * phi(parameters, x, y, z)),
                        )
                    ) - (
                        n2(state, x, y, z) + real(
                            1im *
                            m *
                            bhat(state, parameters, x, y, z) *
                            exp(1im * phi(parameters, x, y, z)),
                        )
                    )
                ),
            ),
        ),
    )

    tracer = TracerNamelist(;
        tracer_setup = :TracerOn,
        leading_order_impact = true,
        initial_tracer = (x, y, z) ->
            n2(state, x, y, z) == 0.0 ? 0.0 :
            real(
                bhat(state, parameters, x, y, z) / n2(state, x, y, z) *
                exp(1im * phi(parameters, x, y, z)),
            ) + z,
        background_tracer = (x, y, z) -> z,
    )

    integrate(
        Namelists(;
            atmosphere,
            domain,
            output,
            tracer,
            # sponge,
            turbulence,
            discretization,
        ),
    )
    return
end