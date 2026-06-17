"""
```julia 
compute_turbulence_diffusion!(state::State)
```

Compute the turbulent eddy diffusivity coefficients by dispatching to a turbulence-scheme-specific method.

```julia
compute_turbulence_diffusion!(state::State, turbulence_scheme::NoTurbulence)
```

Return in case of no turbulence parameterization.

```julia 
compute_turbulence_diffusion!(
    state::State,
    turbulence_scheme::TKEScheme,
)
```

Compute the turbulent eddy diffusivity coefficients for turbulence parameterization using a TKE-Scheme.

The eddy diffusion coefficients for momentum, heat, and turbulence energy are given by 

```math 
\\begin{align*}
    K_M & = L \\sqrt{e_k} \\;, \\\\
    K_H & = L \\sqrt{e_k} \\;, \\\\
    K_{e_k} & = L \\sqrt{e_k} \\;, \\\\
\\end{align*}
```
respectively, with turbulence length scale `L`.
"""
function compute_turbulence_diffusion! end

function compute_turbulence_diffusion!(state::State)
    (; turbulence_scheme) = state.namelists.turbulence

    @dispatch_turbulence_scheme compute_turbulence_diffusion!(
        state, 
        Val(turbulence_scheme),
        )
    return
end

function compute_turbulence_diffusion!(
    state::State,
    turbulence_scheme::Val{:NoTurbulence},
)
    return
end

function compute_turbulence_diffusion!(
    state::State,
    turbulence_scheme::Val{:TKEScheme},
)
    (; model) = state.namelists.atmosphere

    return @dispatch_model compute_turbulence_diffusion!(
        state, turbulence_scheme, Val(model))
end

function compute_turbulence_diffusion!(
    state::State,
    turbulence_scheme::Val{:TKEScheme},
    model::Union{Val{:PseudoIncompressible}, Val{:Compressible}},
)
    (; lv, lb, lt) = state.turbulence.turbulenceconstants
    (; kh, km, kek) = state.turbulence.turbulencediffusioncoefficients
    (; tke) = state.turbulence.turbulencepredictands
    (; rho) = state.variables.predictands
    (; rhobar) = state.atmosphere
    (; k0, k1, j0, j1, i0, i1) = state.domain

    check_tke!(state)
    kh[i0:i1, j0:j1, k0:k1] .=
        lb .*
        sqrt.(
            tke[i0:i1, j0:j1, k0:k1] ./
            (rho[i0:i1, j0:j1, k0:k1] .+ rhobar[i0:i1, j0:j1, k0:k1])
        )
    km[i0:i1, j0:j1, k0:k1] .=
        lv .*
        sqrt.(
            tke[i0:i1, j0:j1, k0:k1] ./
            (rho[i0:i1, j0:j1, k0:k1] .+ rhobar[i0:i1, j0:j1, k0:k1])
        )
    kek[i0:i1, j0:j1, k0:k1] .=
        lt .*
        sqrt.(
            tke[i0:i1, j0:j1, k0:k1] ./
            (rho[i0:i1, j0:j1, k0:k1] .+ rhobar[i0:i1, j0:j1, k0:k1])
        )

    set_boundaries!(state, BoundaryDiffusionCoefficients())
    return
end

function compute_turbulence_diffusion!(
    state::State,
    turbulence_scheme::Val{:TKEScheme},
    model::Val{:Boussinesq},
)
    (; lv, lb, lt) = state.turbulence.turbulenceconstants
    (; kh, km, kek) = state.turbulence.turbulencediffusioncoefficients
    (; tke) = state.turbulence.turbulencepredictands
    (; rhop) = state.variables.predictands
    (; rhobar) = state.atmosphere
    (; k0, k1, j0, j1, i0, i1) = state.domain

    check_tke!(state)
    kh[i0:i1, j0:j1, k0:k1] .=
        lb .*
        sqrt.(
            2.0 .* tke[i0:i1, j0:j1, k0:k1] ./
            (rhop[i0:i1, j0:j1, k0:k1] .+ rhobar[i0:i1, j0:j1, k0:k1])
        )
    km[i0:i1, j0:j1, k0:k1] .=
        lv .*
        sqrt.(
            2.0 .* tke[i0:i1, j0:j1, k0:k1] ./
            (rhop[i0:i1, j0:j1, k0:k1] .+ rhobar[i0:i1, j0:j1, k0:k1])
        )
    kek[i0:i1, j0:j1, k0:k1] .=
        lt .*
        sqrt.(
            2.0 .* tke[i0:i1, j0:j1, k0:k1] ./
            (rhop[i0:i1, j0:j1, k0:k1] .+ rhobar[i0:i1, j0:j1, k0:k1])
        )

    set_boundaries!(state, BoundaryDiffusionCoefficients())
    return
end