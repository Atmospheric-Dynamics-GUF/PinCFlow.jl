# Research/wkb_wp_3d.jl

function wp(;
    x_size::Integer = 128,
    y_size::Integer = 32,
    z_size::Integer = 2400,
    npx::Integer = 1,
    npy::Integer = 1,
    npz::Integer = 1,
    lx::AbstractFloat = 900.0e3,
    ly::AbstractFloat = 30.0e3,
    lz::AbstractFloat = 80.0e3,
    k::AbstractFloat = 0.0,
    l::AbstractFloat = 2 * pi / 30.0e3,
    m::AbstractFloat = 2 * pi / 1.0e3,
    rx::AbstractFloat = 150.0e3,
    ry::AbstractFloat = 0.0,
    rz::AbstractFloat = 5.0e3,
    x0::AbstractFloat = 0.0,
    y0::AbstractFloat = 0.0,
    z0::AbstractFloat = 30.0e3,
    a0::AbstractFloat = 0.7,
    branch::Integer = 1,
    version::AbstractFloat = 2.0,
    coriolis_frequency::AbstractFloat = 2.0e-3,
    output_file::AbstractString = "inertia_gw_wr.h5",
    output_interval::AbstractFloat = 2000.0,
    tmax::AbstractFloat = 2000.0,
    model::Symbol = :Compressible,
    background::Symbol = :Isothermal,
    turbulence_scheme::Symbol = :TKEScheme,
    momentum_coupling::Bool = true,
    entropy_coupling::Bool = true,
    turbulence_filter_order::Integer = 3,
    tracer_setup::Symbol = :TracerOn,
    prepare_restart::Bool = false,
    restart::Bool = false,
    input_file::AbstractString = "inertia_gw_wr_restart.h5",
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
        branch = branch,
    )

    atmosphere = AtmosphereNamelist(; background, model, coriolis_frequency)

    domain = DomainNamelist(;
        x_size,
        y_size,
        z_size,
        lx,
        ly,
        lz,
        npx,
        npy,
        npz,
        nbx = max(3, turbulence_filter_order),
        nby = max(3, turbulence_filter_order),
        nbz = max(3, turbulence_filter_order),
    )

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
        prepare_restart,
        restart,
        iin = -1,
        input_file,
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

    turbulence = TurbulenceNamelist(;
        turbulence_scheme,
        momentum_coupling,
        entropy_coupling,
        turbulence_filter_order,
        turbulence_filter_type = :BoxFilter,
        smooth_turbulence = false,
    )

    discretization = DiscretizationNamelist(; dtmax = 100.0)

    function chi_largescale(z)
        return exp(-(z-z0)^2 / (2 * 5.0e3^2))
    end

    tracer = TracerNamelist(;
        tracer_setup,
        initial_chi = (x, y, z) ->
            -real(
                1im / omega(state, parameters, x, y, z) *
                what(state, parameters, x, y, z) *
                (-(z-z0)/5.0e3^2 * chi_largescale(z)) *
                exp(1im * phi(parameters, x, y, z)),
            ) + chi_largescale(z),
    )

    integrate(
        Namelists(;
            atmosphere,
            domain,
            output,
            discretization,
            turbulence,
            tracer,
        ),
    )

    return
end
