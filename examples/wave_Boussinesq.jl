using Pkg

Pkg.activate("examples")

using PinCFlow

function wave_Boussinesq(;
    output_file::AbstractString = "wave_Boussinesq.h5",
    prepare_restart::Bool = false,
    npx::Integer = 1,
    npy::Integer = 1,
    npz::Integer = 1,
    x_size::Integer = 256,
    y_size::Integer = 256,
    z_size::Integer = 64,
    vertical_boundary_condition::Symbol = :Periodic,
)
    lx = 30.0e3
    ly = 30.0e3
    lz = 1.0e3

    parameters = (
        k = 2 * pi / lx,
        l = 0.0,
        m = 2 * pi / lz,
        rx = 0.0,
        ry = 0.0,
        rz = 0.0, # 5.0e3,
        x0 = 0.0,
        y0 = 0.0,
        z0 = 15.0e3,
        a0 = 0.7,
        version = 2,
    )
    (; k, l, m, a0) = parameters

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
        output_interval = 100.0,
        output_variables = [:thetap, :w, :u, :aux],
        prepare_restart,
        tmax = 2000.0,
    )

    poisson = PoissonNamelist(; initial_cleaning = true)

    turbulence = TurbulenceNamelist(; turbulence_scheme = :NoTurbulence)

    tracer = TracerNamelist(; tracer_setup = :NoTracer)

    integrate(
        Namelists(;
            atmosphere,
            discretization,
            turbulence,
            domain,
            output,
            tracer,
            poisson,
        ),
    )

    return
end
