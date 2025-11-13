"""
```julia
State{
    A <: Namelists,
    B <: Time,
    C <: Constants,
    D <: Domain,
    E <: Grid,
    F <: Atmosphere,
    G <: Sponge,
    H <: Poisson,
    I <: Variables,
    J <: WKB,
    K <: Tracer,
    L <: Ice,
}
```

Model state container.

An instance of this composite type holds complete information about the model configuration and simulation state, so that it is sufficient as primary input to most methods. The construction of such an instance is the first operation performed in [`PinCFlow.Integration.integrate`](@ref), since it almost fully initializes the model.

```julia
State(namelists::Namelists)::State
```

Construct a `State` instance and thus initialize the model.

This method first uses the parameters specified in `namelists` to construct instances of the composite types defined in `FoundationalTypes` (i.e. `Constants`, `Time`, `Domain`, `Grid`, `Atmosphere` and `Sponge`). It then uses these instances to prepare the arrays needed for the Poisson solver, the time integration and the parameterization of unresolved gravity waves with MSGWaM. Afterwards, only three operations of the initialization process remain (these are performed by [`PinCFlow.Integration.integrate`](@ref)), namely the initial cleaning, the setting of the initial ray-volume properties and the reading of input data in restart simulations.

# Fields

  - `namelists::A`: Namelists with all model parameters.

  - `time::B`: Runge-Kutta time integration coefficients.

  - `constants::C`: Physical constants and reference values.

  - `domain::D`: Collection of domain-decomposition and MPI-communication parameters.

  - `grid::E`: Collection of parameters and fields that describe the grid.

  - `atmosphere::F`: Atmospheric-background fields.

  - `sponge::G`: Sponge parameters and damping coefficients.

  - `poisson::H`: Workspace and solution arrays for the Poisson solver.

  - `variables::I`: Arrays needed for the predictions of the prognostic variables.

  - `wkb::J`: Container for WKB ray-tracing data and parameters.

  - `tracer::K`: Tracer setup and parameters.

# Arguments

  - `namelists`: Namelists with all model parameters.

# See also

  - [`PinCFlow.Types.FoundationalTypes.Constants`](@ref)

  - [`PinCFlow.Types.FoundationalTypes.Time`](@ref)

  - [`PinCFlow.Types.FoundationalTypes.Domain`](@ref)

  - [`PinCFlow.Types.FoundationalTypes.Grid`](@ref)

  - [`PinCFlow.Types.FoundationalTypes.Atmosphere`](@ref)

  - [`PinCFlow.Types.FoundationalTypes.Sponge`](@ref)

  - [`PinCFlow.Types.PoissonTypes.Poisson`](@ref)

  - [`PinCFlow.Types.VariableTypes.Variables`](@ref)

  - [`PinCFlow.Types.WKBTypes.WKB`](@ref)

  - [`PinCFlow.Types.TracerTypes.Tracer`](@ref)
"""
struct State{
    A <: Namelists,
    B <: Time,
    C <: Constants,
    D <: Domain,
    E <: Grid,
    F <: Atmosphere,
    G <: Sponge,
    H <: Poisson,
    I <: Variables,
    J <: WKB,
    K <: Tracer,
    L <: Ice,
}
    namelists::A
    time::B
    constants::C
    domain::D
    grid::E
    atmosphere::F
    sponge::G
    poisson::H
    variables::I
    wkb::J
    tracer::K
    ice::L
end

function State(namelists::Namelists)::State

    # Initialize everything.
    constants = Constants(namelists)
    time = Time()
    domain = Domain(namelists)
    grid = Grid(namelists, constants, domain)
    atmosphere = Atmosphere(namelists, constants, domain, grid)
    sponge = Sponge(namelists, domain, grid)
    poisson = Poisson(domain)
    variables = Variables(namelists, constants, domain, atmosphere, grid)
    wkb = WKB(namelists, constants, domain, grid)
    tracer = Tracer(namelists, constants, domain, atmosphere, grid, variables)
    ice = Ice(namelists, constants, domain, atmosphere, grid, variables)

    return State(
        namelists,
        time,
        constants,
        domain,
        grid,
        atmosphere,
        sponge,
        poisson,
        variables,
        wkb,
        tracer,
        ice,
    )
end
