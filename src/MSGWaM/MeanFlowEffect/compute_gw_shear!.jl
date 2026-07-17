"""
```julia
compute_gw_shear!(
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

Compute the gravity-wave induced shear by dispatching to the appropriate method.

```julia
compute_gw_shear!(
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
compute_gw_shear!(
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

Compute the gravity-wave induced shear at ``(i, j, k)``, using

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

  - `wadr`: Phase-space wave-action density.

  - `i`: Zonal grid-cell index.

  - `j`: Meridional grid-cell index.

  - `k`: Vertical grid-cell index.

"""
function compute_gw_shear! end

function compute_gw_shear!(
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

    @dispatch_turbulence_scheme compute_gw_shear!(
        state::State,
        Val(turbulence_scheme),
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

function compute_gw_shear!(
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

@ivy function compute_gw_shear!(
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
    (; rhobar, thetabar) = state.atmosphere
    (; gw_shear) = state.turbulence.turbulencewkbtendencies

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