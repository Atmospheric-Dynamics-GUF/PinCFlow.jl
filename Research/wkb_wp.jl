# Research/wkb_wp_3d.jl

function wkb_wp(;
    x_size::Integer = 32,
    y_size::Integer = 1,
    z_size::Integer = 100,
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
    output_file::AbstractString = "inertia_gw_wkb.h5",
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
    leading_order_impact::Bool = true,
    next_order_impact::Bool = true,
    turbulence_impact::Bool = true,
    restart::Bool = false,
    input_file::AbstractString = "inertia_gw_wkb_restart.h5",
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

    domain = DomainNamelist(; x_size, y_size, z_size, lx, ly, lz, npx, npy, npz)

    atmosphere = AtmosphereNamelist(; background, model, coriolis_frequency)

    state = State(Namelists(; atmosphere, domain))

    atmosphere = AtmosphereNamelist(;
        background,
        model,
        coriolis_frequency,
        initial_u = (
            x,
            y,
            z,
        )->k * wave_action_density(state, parameters, x, y, z) /
           rhobar(state, x, y, z),
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
            :dtkedt,
            :dchidt0,
            :dchidt1,
            :uchi0,
            :vchi0,
            :wchi0,
            :qchi,
            :dchidtq,
            :shear_production,
            :buoyancy_production,
            :gwshear,
            :e,
        ],
        output_interval,
        tmax,
    )

    wkb = WKBNamelist(;
        wkb_mode = :MultiColumn,
        use_saturation = false,
        branch,
        filter_type = :BoxFilter,
        filter_order = 1,
        #nrz = 40,
        #m_bins = 41,
        initial_wave_field = (alpha, x, y, z) -> (
            k,
            l,
            m,
            omega(state, parameters, x, y, z),
            wave_action_density(state, parameters, x, y, z),
        ),
        turbulent_damping = true,
    )

    turbulence = TurbulenceNamelist(;
        turbulence_scheme,
        momentum_coupling,
        entropy_coupling,
        turbulence_filter_order,
        turbulence_filter_type = :BoxFilter,
        smooth_turbulence = false,
    )

    tracer = TracerNamelist(;
        tracer_setup,
        initial_chi = (x, y, z) -> exp(-(z-z0)^2 / (2 * 5.0e3^2)),
        leading_order_impact,
        next_order_impact,
        turbulence_impact,
        apply_lhs_sponge_to_tracer = false,
    )

    discretization = DiscretizationNamelist(; dtmax = 100)

    integrate(
        Namelists(;
            atmosphere,
            domain,
            output,
            wkb,
            discretization,
            turbulence,
            tracer,
        ),
    )

    return
end
