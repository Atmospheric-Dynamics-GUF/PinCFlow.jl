# src/Examples/wkb_boussinesq_periodic.jl

function wkb_boussinesq_periodic(;
    output_file::AbstractString = "wkb_boussinesq_periodic.h5",
    x_size::Integer = 31,
    y_size::Integer = 1,
    z_size::Integer = 31,
    lambdah::AbstractFloat = 300.0e3,
    a0::AbstractFloat = 1.0,
)
    lx = 4000.0e3
    ly = 300.0e3
    lz = 20.0e3

    parameters = (
        k = 0.0,
        l = 2 * pi / lambdah,
        m = 2 * pi / 1.0e3,
        a0 = a0,
        lx = lx,
        ly = ly,
        lz = lz,
    )
    (; k, l, m) = parameters

    model = :Boussinesq
    background = :StableStratification
    coriolis_frequency = 1.0e-4
    buoyancy_frequency = 2.0e-2

    domain = DomainNamelist(;
        lx,
        ly,
        lz,
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

    wkb = WKBNamelist(;
        wkb_mode = :MultiColumn,
        initial_wave_field = (alpha, x, y, z) -> (
            k,
            l,
            m,
            omega(state, parameters, x, y, z),
            wave_action_density(state, parameters, x, y, z),
        ),
        turbulent_damping = false,
        use_saturation = false,
    )

    output = OutputNamelist(;
        output_file,
        output_interval = 1.0,
        output_variables = [
            :dchidt0,
            :dchidtq,
            :tke,
            :u,
            :e,
            :dtkedt,
            :uchi0,
            :vchi0,
            :wchi0,
            :qchi,
            :tracerdiffusion,
        ],
        tmax = 1.0,
    )

    discretization = DiscretizationNamelist(; dtmax = 100.0)

    function chi_largescale(x::Real, z::Real)
        return (1 + cos(2 * pi * (z - lz / 2) / lz)) / 2
    end

    tracer = TracerNamelist(;
        tracer_setup = :TracerOn,
        initial_chi = (x, y, z) -> chi_largescale(x, z),
        leading_order_impact = true,
        turbulence_wave_impact = true,
    )

    lt = 30.0

    turbulence = TurbulenceNamelist(;
        gw_coupling = true,
        initial_tke = (x, y, z) ->
            max(
                0,
                lt * (
                    lt * m^2 / 2 * (
                        abs(uhat(state, parameters, x, y, z))^2 +
                        abs(vhat(state, parameters, x, y, z))^2
                    ) - lt * (n2(state, x, y, z))
                ),
            ) / 2,
        momentum_coupling = false,
    )

    integrate(
        Namelists(;
            atmosphere,
            wkb,
            domain,
            output,
            turbulence,
            tracer,
            discretization,
        ),
    )

    return
end
