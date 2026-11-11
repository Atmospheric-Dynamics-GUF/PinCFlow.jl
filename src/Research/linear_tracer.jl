# src/Research/linear_tracer.jl

function linear_tracer(;
    npx::Integer = 1,
    npy::Integer = 1,
    npz::Integer = 1,
    output_file::AbstractString = "linear_tracer.h5",
    x_size::Integer = 214,
    y_size::Integer = 16,
    z_size::Integer = 200,
    a0::AbstractFloat = 1.0,
)
    lx = 4000.0e3
    ly = 300.0e3
    lz = 30.0e3

    parameters = (
        k = 0.0,
        l = 2 * pi / ly,
        m = 2 * pi / 1.0e3,
        a0 = a0,
        x0 = 0.0,
        y0 = 0.0,
        z0 = 1000.0e3,
        rx = lx / 2,
        ry = 0.0,
        rz = 5.0e3,
    )

    model = :Boussinesq
    background = :StableStratification
    coriolis_frequency = 1.0e-4
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
        output_interval = 3600.0,
        output_variables = [:u, :v, :w, :rhop],
        tmax = 3600.0 * 25,
    )

    alpha = 1.0

    tracer = TracerNamelist(;
        tracer_setup = :TracerOn,
        background_tracer_gradient = alpha,
        initial_chi = (x, y, z) -> real(
            -1im / omega(state, parameters, x, y, z) *
            alpha *
            what(state, parameters, x, y, z),
        ),
    )

    turbulence = TurbulenceNamelist(; turbulence_scheme = :NoTurbulence)

    integrate(Namelists(; atmosphere, domain, output, tracer, turbulence))

    return
end
