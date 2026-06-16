# Turbulence_research/wp_1d.jl

function wp_1d(;
    x_size::Integer = 32,
    y_size::Integer = 1,
    z_size::Integer = 960,
    npx::Integer = 1,
    npy::Integer = 1,
    npz::Integer = 1,
    output_file::AbstractString = "wp_1d.h5",
    output_interval::AbstractFloat = 360.0,
    tmax::AbstractFloat = 3600.0,
    a0::AbstractFloat = 1.0,
    model::Symbol = :Compressible,
    background::Symbol = :Isothermal,
    momentum_coupling::Bool = true,
    entropy_coupling::Bool = true,
)
    lx = 300.0e3
    ly = 30.0e3
    lz = 60.0e3

    parameters = (
        k = 2 * pi / lx,
        l = 0.0,
        m = 2 * pi / 3.0e3,
        rx = 0.0,
        ry = 0.0,
        rz = 0.05,
        x0 = 0.0,
        y0 = 0.0,
        z0 = 30.e3,
        a0,
    )

    coriolis_frequency = 0.0

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

    integrate(Namelists(; atmosphere, domain, output, turbulence))
    return
end
