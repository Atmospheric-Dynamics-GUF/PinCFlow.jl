function wave_Boussinesq(;
    output_file::AbstractString = "wave_Boussinesq.h5",
    prepare_restart::Bool = false,
    npx::Integer = 1,
    npy::Integer = 1,
    npz::Integer = 1,
    x_size::Integer = 64,
    y_size::Integer = 1,
    z_size::Integer = 64,
    vertical_boundary_condition::Symbol = :Periodic,
    visualize::Bool = false,
    output_interval::AbstractFloat = 100.0,
    tmax::AbstractFloat = 100.0,
)
    lx = 30.0e3
    ly = 30.0e3
    lz = 3.0e3

    parameters = (
        k = 2 * pi / lx,
        l = 0.0,
        m = 6 * pi / lz,
        rx = 0.0,
        ry = 0.0,
        rz = 0.0,
        x0 = 0.0,
        y0 = 0.0,
        z0 = 0.0,
        a0 = 0.1,
        version = 2,
    )
    (; k, l, m, a0) = parameters

    model = :Boussinesq
    background = :StableStratification
    coriolis_frequency = 0.0
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
        vertical_boundary_condition,
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
        kinematic_viscosity = 0.0,
        thermal_conductivity = 0.0,
        kinematic_diffusivity = 0.0,
    )

    discretization = DiscretizationNamelist(; dtmax = 100)

    output = OutputNamelist(;
        output_file,
        output_interval,
        output_variables = [:rhop, :w, :v, :u, :tke],
        prepare_restart,
        tmax,
    )

    turbulence = TurbulenceNamelist(;
        turbulence_scheme = :NoTurbulence,
        momentum_coupling = true,
        tracer_coupling = true,
    )

    kenv = 2 * pi / lx
    lenv = 0.0
    menv = 2 * pi / lz

    function chils(x::Real, y::Real, z::Real)
        return 1/2^3 *
               (1 + cos(kenv * x)) *
               (1 + cos(lenv * y)) *
               (1 + cos(menv * (z - lz / 2)))
    end

    function dxchils(x::Real, y::Real, z::Real)
        return -kenv/2^3 *
               sin(kenv * x) *
               (1 + cos(lenv * y)) *
               (1 + cos(menv * (z - lz / 2)))
    end

    function dychils(x::Real, y::Real, z::Real)
        return -lenv/2^3 *
               (1 + cos(kenv * x)) *
               sin(lenv * y) *
               (1 + cos(menv * (z - lz / 2)))
    end

    function dzchils(x::Real, y::Real, z::Real)
        return -menv/2^3 *
               (1 + cos(kenv * x)) *
               (1 + cos(lenv * y)) *
               sin(menv * z)
    end

    tracer = TracerNamelist(;
        tracer_setup = :NoTracer,
        initial_chi = (x, y, z) ->
            chils(x, y, z) + real(
                -1im/omega(state, parameters, x, y, z) *
                (
                    uhat(state, parameters, x, y, z) * dxchils(x, y, z) +
                    vhat(state, parameters, x, y, z) * dychils(x, y, z) +
                    what(state, parameters, x, y, z) * dzchils(x, y, z)
                ) *
                exp(1im * phi(parameters, x, y, z)),
            ),
    )

    integrate(
        Namelists(;
            atmosphere,
            discretization,
            turbulence,
            domain,
            output,
            tracer,
        ),
    )

    return
end
