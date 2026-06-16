# Turbulence_research/wkb_wp_3d.jl

function wp(;
    x_size::Integer = 512,
    y_size::Integer = 16,
    z_size::Integer = 1000,
    npx::Integer = 1,
    npy::Integer = 1,
    npz::Integer = 1,
    lx::AbstractFloat = 9000.e3,
    ly::AbstractFloat = 300.e3,
    lz::AbstractFloat = 100.e3,
    k::AbstractFloat = 0.0,
    l::AbstractFloat = 2 * pi / 300.e3,
    m::AbstractFloat = 2 * pi / 1.e3,
    rx::AbstractFloat = 1500.e3,
    ry::AbstractFloat = 0.0,
    rz::AbstractFloat = 5.e3,
    x0::AbstractFloat = 0.0,
    y0::AbstractFloat = 0.0,
    z0::AbstractFloat = 30.e3,
    a0::AbstractFloat = 1.0,
    version::AbstractFloat = 2.0,
    coriolis_frequency::AbstractFloat = 1.e-4,
    output_file::AbstractString = "wp_3d.h5",
    output_interval::AbstractFloat = 360.0,
    tmax::AbstractFloat = 3600.0,
    model::Symbol = :Compressible,
    background::Symbol = :Isothermal,
    momentum_coupling::Bool = true,
    entropy_coupling::Bool = true,
    tracer_setup::Symbol = :TracerOn,
)
    parameters = (
        k = k,
        l = l,
        m = m,
        rx = rx,
        ry = ry,
        rz = rz,
        x0 = x0,
        y0 = y0,
        z0 = z0,
        a0 = a0,
        version = version,
    )

    atmosphere = AtmosphereNamelist(; background, model, coriolis_frequency)

    domain = DomainNamelist(; x_size, y_size, z_size, lx, ly, lz, npx, npy, npz)

    state = State(Namelists(; atmosphere, domain))
    (; g) = state.constants

    atmosphere = AtmosphereNamelist(;
        background,
        model,
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
        output_variables = [
            :u,
            :v,
            :w,
            :rhop,
            :tke,
            :shear_production,
            :buoyancy_production,
        ],
        output_interval,
        tmax,
    )

    turbulence = TurbulenceNamelist(; momentum_coupling, entropy_coupling)

    tracer = TracerNamelist(;
        tracer_setup,
        initial_chi = (x, y, z) ->
            -real(
                1im / omega(state, parameters, x, y, z) *
                what(state, parameters, x, y, z) *
                exp(1im * phi(parameters, x, y, z)),
            ) + z,
    )

    integrate(Namelists(; atmosphere, domain, output, turbulence, tracer))

    return
end
