"""
```julia
TracerWKBTendencies{A <: AbstractArray{<:AbstractFloat, 3}}
```

Tracer impact of unresolved gravity waves.

```julia 
TracerWKBTendencies(
    namelists::Namelists,
    domain::Domain,
)::TracerWKBTendencies
```

Construct a `TracerWKBTendencies` instance by dispatching to the appropriate method.

```julia
TracerWKBTendencies(
    namelists::Namelists,
    domain::Domain,
    tracer_setup::Val{:NoTracer},
)::TracerWKBTendencies
```

Construct a `TracerWKBTendencies` instance with zero-size arrays for configurations without tracer transport.

```julia 
TracerWKBTendencies(
    namelists::Namelists,
    domain::Domain,
    tracer_setup::Val{:TracerOn},
)::TracerWKBTendencies
```

Construct a `TracerWKBTendencies` instance by dispatching to the appropriate method.

```julia 
TracerWKBTendencies(
    namelists::Namelists,
    domain::Domain,
    wkb_mode::Val{:NoWKB},
)::TracerWKBTendencies
```

Construct a `TracerWKBTendencies` instance with zero-size arrays for non-WKB configurations.

```julia
TracerWKBTendencies(
    namelists::Namelists,
    domain::Domain,
    wkb_mode::Union{Val{:SteadyState}, Val{:SingleColumn}, Val{:MultiColumn}},
)::TracerWKBTendencies
```

Construct a `TracerWKBTendencies` instance with zero-initialized arrays if `state.namelists.tracer.leading_order_impact == true`, otherwise the arrays are zero-size.

# Fields 

  - `dchidt0::A`: Leading-order tracer impact of unresolved gravity waves.

  - `dchidt1::B`: Next-order tracer impact of unresolved gravity waves.

  - `dchidtq::C`: Tracer impact of turbulence induced by unresolved gravity waves.

# Arguments

  - `namelists`: Namelists with all model parameters.

  - `domain`: Collection of domain-decomposition and MPI-communication parameters.

  - `tracer_setup`: General tracer-transport configuration.

  - `wkb_mode`: Approximations used by MS-GWaM.
"""
struct TracerWKBTendencies{
    A <: AbstractArray{<:AbstractFloat, 3},
    B <: AbstractArray{<:AbstractFloat, 3},
    C <: AbstractArray{<:AbstractFloat, 3},
}
    dchidt0::A
    dchidt1::B
    dchidtq::C
end

function TracerWKBTendencies(
    namelists::Namelists,
    domain::Domain,
)::TracerWKBTendencies
    (; tracer_setup) = namelists.tracer

    @dispatch_tracer_setup return TracerWKBTendencies(
        namelists,
        domain,
        Val(tracer_setup),
    )
end

function TracerWKBTendencies(
    namelists::Namelists,
    domain::Domain,
    tracer_setup::Val{:NoTracer},
)::TracerWKBTendencies
    return TracerWKBTendencies([zeros(0, 0, 0) for i in 1:3]...)
end

function TracerWKBTendencies(
    namelists::Namelists,
    domain::Domain,
    tracer_setup::Val{:TracerOn},
)::TracerWKBTendencies
    (; wkb_mode) = namelists.wkb

    @dispatch_wkb_mode return TracerWKBTendencies(
        namelists,
        domain,
        Val(wkb_mode),
    )
end

function TracerWKBTendencies(
    namelists::Namelists,
    domain::Domain,
    wkb_mode::Val{:NoWKB},
)::TracerWKBTendencies
    return TracerWKBTendencies([zeros(0, 0, 0) for i in 1:3]...)
end

function TracerWKBTendencies(
    namelists::Namelists,
    domain::Domain,
    wkb_mode::Union{Val{:SteadyState}, Val{:SingleColumn}, Val{:MultiColumn}},
)::TracerWKBTendencies
    (; nxx, nyy, nzz) = domain
    (; leading_order_impact, next_order_impact, turbulence_impact) =
        namelists.tracer

    if leading_order_impact
        dchidt0 = zeros(nxx, nyy, nzz)
    else
        dchidt0 = zeros(0, 0, 0)
    end
    if next_order_impact
        dchidt1 = zeros(nxx, nyy, nzz)
    else
        dchidt1 = zeros(0, 0, 0)
    end
    if turbulence_impact
        dchidtq = zeros(nxx, nyy, nzz)
    else
        dchidtq = zeros(0, 0, 0)
    end

    return TracerWKBTendencies(dchidt0, dchidt1, dchidtq)
end
