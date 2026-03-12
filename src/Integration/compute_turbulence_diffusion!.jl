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

Compute the turbulent eddy diffusivity coefficients for turbulence parameterization using a TKE-Scheme by dispatching over a model-specific method.

```julia 
compute_turbulence_diffusion!(
    state::State,
    turbulence_scheme::TKEScheme,
    model::Union{PseudoIncompressible, Compressible},
)
```

Compute the turbulent eddy diffusivity coefficients for turbulence parameterization using a TKE-Scheme in pseudo-incompressible and compressible mode. 

The eddy diffusion coefficients for momentum, heat, and turbulent kinetic energy are given by 

```math 
\\begin{align*}
    K_M & = l_v \\sqrt{2 e_\\mathrm{k}} \\;, \\\\
    K_H & = l_h \\sqrt{2 e_\\mathrm{k}} \\;, \\\\
    K_{e_\\mathrm{k}} & = l_t \\sqrt{2 e_\\mathrm{k}} \\;, \\\\
\\end{align*}
```
respectively, with turbulence length scales `l_v`, `l_h`, and `l_t` and mass-specific turbulent kinetic energy `e_\\mathrm{k}`.

```julia 
compute_turbulence_diffusion!(
    state::State,
    turbulence_scheme::TKEScheme,
    model::Boussinesq,
)
```

Compute the turbulent eddy diffusivity coefficients for turbulence parameterization using a TKE-Scheme in Boussinesq mode.

The eddy diffusion coefficients for momentum, heat, and turbulent kinetic energy are given by 

```math 
\\begin{align*}
    K_M & = l_v \\sqrt{2 e_\\mathrm{k}} \\;, \\\\
    K_H & = l_h \\sqrt{2 e_\\mathrm{k}} \\;, \\\\
    K_{e_\\mathrm{k}} & = l_t \\sqrt{2 e_\\mathrm{k}} \\;, \\\\
\\end{align*}
```
respectively, with turbulence length scales `l_v`, `l_h`, and `l_t` and mass-specific turbulent kinetic energy `e_\\mathrm{k}`.

# Arguements:

  - `state`: Model state. 

  - `turbulence_scheme`: General turbulence-parameterization configuration.

  - `model`: Dynamic equations.
"""
function compute_turbulence_diffusion! end

function compute_turbulence_diffusion!(state::State)
    (; turbulence_scheme) = state.namelists.turbulence

    compute_turbulence_diffusion!(state, turbulence_scheme)
    return
end

function compute_turbulence_diffusion!(
    state::State,
    turbulence_scheme::NoTurbulence,
)
    return
end

function compute_turbulence_diffusion!(
    state::State,
    turbulence_scheme::TKEScheme,
)
    (; model) = state.namelists.atmosphere

    return compute_turbulence_diffusion!(state, turbulence_scheme, model)
end

function compute_turbulence_diffusion!(
    state::State,
    turbulence_scheme::TKEScheme,
    model::Union{PseudoIncompressible, Compressible},
)
    (; lv, lb, lt, prandtlinv) = state.turbulence.turbulenceconstants
    (; kh, km, kek) = state.turbulence.turbulencediffusioncoefficients
    (; tke) = state.turbulence.turbulencepredictands
    (; rho) = state.variables.predictands
    (; rhobar) = state.atmosphere
    (; k0, k1, j0, j1, i0, i1) = state.domain

    check_tke!(state)
    kh[i0:i1, j0:j1, k0:k1] .=
        lb .*
        sqrt.(
            2.0 * tke[i0:i1, j0:j1, k0:k1] ./
            (rho[i0:i1, j0:j1, k0:k1] .+ rhobar[i0:i1, j0:j1, k0:k1])
        )
    km[i0:i1, j0:j1, k0:k1] .=
        lv .*
        sqrt.(
            2.0 * tke[i0:i1, j0:j1, k0:k1] ./
            (rho[i0:i1, j0:j1, k0:k1] .+ rhobar[i0:i1, j0:j1, k0:k1])
        )
    kek[i0:i1, j0:j1, k0:k1] .=
        lt .*
        sqrt.(
            2.0 * tke[i0:i1, j0:j1, k0:k1] ./
            (rho[i0:i1, j0:j1, k0:k1] .+ rhobar[i0:i1, j0:j1, k0:k1])
        )

    set_boundaries!(state, BoundaryDiffusionCoefficients())
    return
end

function compute_turbulence_diffusion!(
    state::State,
    turbulence_scheme::TKEScheme,
    model::Boussinesq,
)
    (; lb, lv, lt, prandtlinv) = state.turbulence.turbulenceconstants
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
