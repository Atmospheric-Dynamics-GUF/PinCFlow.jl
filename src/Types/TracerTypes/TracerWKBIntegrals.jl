"""
```julia
TracerWKBIntegrals{A <: AbstractArray{<:AbstractFloat, 3}}
```

Integrals of gravity-wave induced tracer fluxes.

```julia 
TracerWKBIntegrals(
    namelists::Namelists,
    domain::Domain,
)::TracerWKBIntegrals
```

Construct a `TracerWKBIntegrals` instance by dispatching to the appropriate method.

```julia
TracerWKBIntegrals(
    namelists::Namelists,
    domain::Domain,
    tracer_setup::Val{:NoTracer},
)::TracerWKBIntegrals
```

Construct a `TracerWKBIntegrals` instance with zero-size arrays for configurations without tracer transport.

```julia 
TracerWKBIntegrals(
    namelists::Namelists,
    domain::Domain,
    tracer_setup::Val{:TracerOn},
)::TracerWKBIntegrals
```

Construct a `TracerWKBIntegrals` instance by dispatching to the appropriate method.

```julia 
TracerWKBIntegrals(
    namelists::Namelists,
    domain::Domain,
    wkb_mode::Val{:NoWKB},
)::TracerWKBIntegrals
```

Construct a `TracerWKBIntegrals` instance with zero-size arrays for non-WKB configurations.

```julia
TracerWKBIntegrals(
    namelists::Namelists,
    domain::Domain,
    wkb_mode::Union{Val{:SteadyState}, Val{:SingleColumn}, Val{:MultiColumn}},
)::TracerWKBIntegrals
```

Construct a `TracerWKBIntegrals` instance with zero-initialized arrays if `state.namelists.tracer.leading_order_impact == true`, otherwise the arrays are zero-size.

# Fields 

  - `uchi0::A`: Leading-order zonal tracer fluxes.

  - `vchi0::A`: Leading-order meridional tracer fluxes.

  - `wchi0::A`: Leading-order vertical tracer fluxes.

  - `uchi1::A`: Next-order zonal tracer fluxes.

  - `vchi1::A`: Next-order meridional tracer fluxes.

  - `wchi1::A`: Next-order vertical tracer fluxes.

  - `uhat::B`: Leading-order zonal wind amplitude of unresolved gravity waves.

  - `vhat::B`: Leading-order meridional wind amplitude of unresolved gravity waves.

  - `what::B`: Leading-order vertical wind amplitude of unresolved gravity waves.

  - `bhat::B`: Leading-order buoyancy amplitude of unresolved gravity waves.

  - `pihat::B`: Leading-order Exner pressure amplitude of unresolved gravity waves.

  - `chihat::B`: Leading-order tracer amplitude of unresolved gravity waves.

  - `uhatold::B`: Leading-order zonal wind amplitude of unresolved gravity waves of the previous time step.

  - `vhatold::B`: Leading-order meridional wind amplitude of unresolved gravity waves of the previous time step.

  - `whatold::B`: Leading-order vertical wind amplitude of unresolved gravity waves of the previous time step.

  - `bhatold::B`: Leading-order buoyancy amplitude of unresolved gravity waves of the previous time step.

  - `chihatold::B`: Leading-order tracer amplitude of unresolved gravity waves of the previous time step.

# Arguments

  - `namelists`: Namelists with all model parameters.

  - `domain`: Collection of domain-decomposition and MPI-communication parameters.

  - `tracer_setup`: General tracer-transport configuration.

  - `wkb_mode`: Approximations used by MS-GWaM.
"""
struct TracerWKBIntegrals{
    A <: AbstractArray{<:AbstractFloat, 3},
    B <: AbstractArray{<:ComplexF64, 3},
}
    uchi0::A
    vchi0::A
    wchi0::A
    uchi1::A
    vchi1::A
    wchi1::A
    uhat::B
    vhat::B
    what::B
    bhat::B
    pihat::B
    chihat::B
    uhatold::B
    vhatold::B
    whatold::B
    bhatold::B
    chihatold::B
end

function TracerWKBIntegrals(
    namelists::Namelists,
    domain::Domain,
)::TracerWKBIntegrals
    (; tracer_setup) = namelists.tracer

    @dispatch_tracer_setup return TracerWKBIntegrals(
        namelists,
        domain,
        Val(tracer_setup),
    )
end

function TracerWKBIntegrals(
    namelists::Namelists,
    domain::Domain,
    tracer_setup::Val{:NoTracer},
)::TracerWKBIntegrals
    return TracerWKBIntegrals(
        [zeros(0, 0, 0) for i in 1:6]...,
        [zeros(ComplexF64, 0, 0, 0) for i in 1:11]...,
    )
end

function TracerWKBIntegrals(
    namelists::Namelists,
    domain::Domain,
    tracer_setup::Val{:TracerOn},
)::TracerWKBIntegrals
    (; wkb_mode) = namelists.wkb

    @dispatch_wkb_mode return TracerWKBIntegrals(
        namelists,
        domain,
        Val(wkb_mode),
    )
end

function TracerWKBIntegrals(
    namelists::Namelists,
    domain::Domain,
    wkb_mode::Val{:NoWKB},
)::TracerWKBIntegrals
    return TracerWKBIntegrals(
        [zeros(0, 0, 0) for i in 1:6]...,
        [zeros(ComplexF64, 0, 0, 0) for i in 1:11]...,
    )
end

function TracerWKBIntegrals(
    namelists::Namelists,
    domain::Domain,
    wkb_mode::Union{Val{:SteadyState}, Val{:SingleColumn}, Val{:MultiColumn}},
)::TracerWKBIntegrals
    (; nxx, nyy, nzz) = domain
    (; leading_order_impact, next_order_impact) = namelists.tracer

    if leading_order_impact
        nxl = nxx
        nyl = nyy
        nzl = nzz
    else
        nxl = 0
        nyl = 0
        nzl = 0
    end

    if next_order_impact
        nxn = nxx
        nyn = nyy
        nzn = nzz
    else
        nxn = 0
        nyn = 0
        nzn = 0
    end

    return TracerWKBIntegrals(
        [zeros(nxl, nyl, nzl) for i in 1:3]...,
        [zeros(nxn, nyn, nzn) for i in 1:3]...,
        [zeros(ComplexF64, nxn, nyn, nzn) for i in 1:11]...,
    )
end
