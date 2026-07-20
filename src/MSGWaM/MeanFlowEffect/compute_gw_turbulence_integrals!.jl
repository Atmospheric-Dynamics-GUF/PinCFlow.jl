"""
```julia
compute_gw_turbulence_integrals!(
    state::State,
    fc::AbstractFloat,
    omir::AbstractFloat,
    kr::AbstractFloat,
    lr::AbstractFloat,
    mr::AbstractFloat,
    wadr::AbstractFloat,
    i::Integer,
    j::Integer,
    k::Integer,
)
```

Compute the gravity-wave shear by dispatching to the appropriate method.

```julia
compute_gw_turbulence_integrals!(
    state::State,
    turbulence_scheme::Val{:NoTurbulence},
    fc::AbstractFloat,
    omir::AbstractFloat,
    kr::AbstractFloat,
    lr::AbstractFloat,
    mr::AbstractFloat,
    wadr::AbstractFloat,
    i::Integer,
    j::Integer,
    k::Integer,
)
```

Return for configurations without turbulence parameterization.

```julia
compute_gw_turbulence_integrals!(
    state::State,
    turbulence_scheme::Val{:TKEScheme},
    fc::AbstractFloat,
    omir::AbstractFloat,
    kr::AbstractFloat,
    lr::AbstractFloat,
    mr::AbstractFloat,
    wadr::AbstractFloat,
    i::Integer,
    j::Integer,
    k::Integer,
)
```

Compute the gravity-wave shear at ``(i, j, k)``, using

```math
\\begin{align*}
    \\mathcal{S}_{gw} = \\frac{m_r^4}{\\bar{\\rho}} \\frac{\\left(f^2 + \\hat{\\omega}_r^2\\right)}{\\hat{\\omega}_r \\left|\\boldsymbol{k}_r\\right|^2}.
\\end{align*}
```

# Arguments

  - `state`: Model state.

  - `turbulence_scheme`: General turbulence parameterization configuration.

  - `fc`: Coriolis parameter.

  - `omir`: Gravity-wave intrinsic frequency.

  - `kr`: Zonal wavenumber.

  - `lr`: Meridional wavenumber.

  - `mr`: Vertical wavenumber.

  - `wadr`: Physical-space wave-action density.

  - `i`: Zonal grid-cell index.

  - `j`: Meridional grid-cell index.

  - `k`: Vertical grid-cell index.

!!! danger "Experimental"
    The gravity-wave shear is an experimental feature that hasn't been validated yet.
"""
function compute_gw_turbulence_integrals! end

function compute_gw_turbulence_integrals!(
    state::State,
    fc::AbstractFloat,
    omir::AbstractFloat,
    kr::AbstractFloat,
    lr::AbstractFloat,
    mr::AbstractFloat,
    wadr::AbstractFloat,
    i::Integer,
    j::Integer,
    k::Integer,
)
    (; turbulence_scheme) = state.namelists.turbulence

    @dispatch_turbulence_scheme compute_gw_turbulence_integrals!(
        state,
        Val(turbulence_scheme),
        fc,
        omir,
        kr,
        lr,
        mr,
        wadr,
        i,
        j,
        k,
    )

    return
end

function compute_gw_turbulence_integrals!(
    state::State,
    turbulence_scheme::Val{:NoTurbulence},
    fc::AbstractFloat,
    omir::AbstractFloat,
    kr::AbstractFloat,
    lr::AbstractFloat,
    mr::AbstractFloat,
    wadr::AbstractFloat,
    i::Integer,
    j::Integer,
    k::Integer,
)
    return
end

@ivy function compute_gw_turbulence_integrals!(
    state::State,
    turbulence_scheme::Val{:TKEScheme},
    fc::AbstractFloat,
    omir::AbstractFloat,
    kr::AbstractFloat,
    lr::AbstractFloat,
    mr::AbstractFloat,
    wadr::AbstractFloat,
    i::Integer,
    j::Integer,
    k::Integer,
)
    (; rhobar) = state.atmosphere
    (; gw_shear) = state.turbulence.turbulencewkbintegrals

    gw_shear[i, j, k] +=
        mr^4  * 
        wadr / 
        rhobar[i, j, k] * 
        (fc^2 + omir^2) / (
            omir * 
            (kr^2 + lr^2 + mr^2)
        )
        
    return
end
