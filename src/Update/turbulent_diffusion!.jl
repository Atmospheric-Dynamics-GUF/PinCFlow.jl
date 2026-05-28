"""
```julia 
turbulent_diffusion!(state::State, dt::AbstractFloat)
```

Apply diffusion to the momentum, mass-weighted potential temperature, and tracers by dispatching to turbulence parameterization-specific method.

```julia 
turbulent_diffusion!(
    state::State,
    dt::AbstractFloat,
    turbulence_scheme::Val{:NoTurbulence},
)
```

Return for configurations without turbulence parameterization.

```julia 
turbulent_diffusion!(
    state::State,
    dt::AbstractFloat,
    turbulence_scheme::Val{:TKEScheme},
)
```

Apply diffusion by dispatching to specialized methods for momentum, mass-weighted potential temperature, and tracers, according to configurations set by `state.namelist.turbulence.momentum_coupling`, `state.namelist.turbulence.entropy_coupling`, and `state.namelist.turbulence.tracer_coupling`, respectively.

```julia 
turbulent_diffusion!(state::State, dt::AbstractFloat, variable::U)
```

Apply diffusion to the zonal momentum. 

The prognostic equation

```math
\\frac{\\partial u}{\\partial t} = \\frac{1}{J}\\frac{\\partial}{\\partial \\hat{z}}\\left(\\frac{K_\\mathrm{M}}{J}\\frac{\\partial u}{\\partial \\hat{z}}\\right)
```

is solved using the Crank-Nicolson scheme, where the system 

```math
a_{i+1/2,k} u_{i+1/2,k-1}^{n+1} + b_{i+1/2,k} u_{i+1/2,k}^{n+1} + c_{i+1/2,k} u_{i+1/2,k+1}^{n+1} = f_{i+1/2,k}
```

is solved using a Thomas tridiagonal solver, with ``\\mathcal{K}_\\mathrm{M} = \\frac{K_\\mathrm{M}}{J}`` and

```math 
\\begin{align*}
    a_{i+1/2,k} &= -\\frac{\\Delta t}{2(\\Delta \\hat{z})^2}\\frac{\\mathcal{K}_{\\mathrm{M},i+1/2,k-1/2}}{J_{i+1/2}} ,\\\\
    b_{i+1/2,k} &= 1 + \\frac{\\Delta t}{2(\\Delta \\hat{z})^2}\\frac{\\mathcal{K}_{\\mathrm{M},i+1/2,k+1/2}}{J_{i+1/2}} 
        + \\frac{\\Delta t}{2(\\Delta \\hat{z})^2}\\frac{\\mathcal{K}_{\\mathrm{M},i+1/2,k-1/2}}{J_{i+1/2}} ,\\\\
    c_{i+1/2,k} &= -\\frac{\\Delta t}{2(\\Delta \\hat{z})^2}\\frac{\\mathcal{K}_{\\mathrm{M},i+1/2,k+1/2}}{J_{i+1/2}} ,\\\\
    f_{i+1/2,k} &= \\left[1 - \\frac{\\Delta t}{2(\\Delta \\hat{z})^2}\\frac{\\mathcal{K}_{\\mathrm{M},i+1/2,k+1/2}}{J_{i+1/2}} 
        - \\frac{\\Delta t}{2(\\Delta \\hat{z})^2}\\frac{\\mathcal{K}_{\\mathrm{M},i+1/2,k-1/2}}{J_{i+1/2}}\\right]v_{i+1/2}^{n} \\\\
        &+\\frac{\\Delta t}{2(\\Delta \\hat{z})^2}\\frac{\\mathcal{K}_{\\mathrm{M},i+1/2,k+1/2}}{J_{i+1/2}}v_{i+1/2,k+1}^{n}
        +\\frac{\\Delta t}{2(\\Delta \\hat{z})^2}\\frac{\\mathcal{K}_{\\mathrm{M},i+1/2,k-1/2}}{J_{i+1/2}}v_{i+1/2,k-1}^{n} .
\\end{align*}
```

```julia 
turbulent_diffusion!(state::State, dt::AbstractFloat, variable::V)
```

Apply diffusion to the meridional momentum. 

The prognostic equation

```math
\\frac{\\partial v}{\\partial t} = \\frac{1}{J}\\frac{\\partial}{\\partial \\hat{z}}\\left(\\frac{K_\\mathrm{M}}{J}\\frac{\\partial v}{\\partial \\hat{z}}\\right)
```

is solved using the Crank-Nicolson scheme, where the system 

```math
a_{j+1/2,k} v_{j+1/2,k-1}^{n+1} + b_{j+1/2,k} v_{j+1/2,k}^{n+1} + c_{j+1/2,k} v_{j+1/2,k+1}^{n+1} = f_{j+1/2,k}
```

is solved using a Thomas tridiagonal solver, with ``\\mathcal{K}_\\mathrm{M} = \\frac{K_\\mathrm{M}}{J}`` and

```math 
\\begin{align*}
    a_{j+1/2,k} &= -\\frac{\\Delta t}{2(\\Delta \\hat{z})^2}\\frac{\\mathcal{K}_{\\mathrm{M},j+1/2,k-1/2}}{J_{j+1/2}} ,\\\\
    b_{j+1/2,k} &= 1 + \\frac{\\Delta t}{2(\\Delta \\hat{z})^2}\\frac{\\mathcal{K}_{\\mathrm{M},j+1/2,k+1/2}}{J_{j+1/2}} 
        + \\frac{\\Delta t}{2(\\Delta \\hat{z})^2}\\frac{\\mathcal{K}_{\\mathrm{M},j+1/2,k-1/2}}{J_{j+1/2}} ,\\\\
    c_{j+1/2,k} &= -\\frac{\\Delta t}{2(\\Delta \\hat{z})^2}\\frac{\\mathcal{K}_{\\mathrm{M},j+1/2,k+1/2}}{J_{j+1/2}} ,\\\\
    f_{j+1/2,k} &= \\left[1 - \\frac{\\Delta t}{2(\\Delta \\hat{z})^2}\\frac{\\mathcal{K}_{\\mathrm{M},j+1/2,k+1/2}}{J_{j+1/2}} 
        - \\frac{\\Delta t}{2(\\Delta \\hat{z})^2}\\frac{\\mathcal{K}_{\\mathrm{M},j+1/2,k-1/2}}{J_{j+1/2}}\\right]v_{j+1/2}^{n}\\\\
        &+\\frac{\\Delta t}{2(\\Delta \\hat{z})^2}\\frac{\\mathcal{K}_{\\mathrm{M},j+1/2,k+1/2}}{J_{j+1/2}}v_{j+1/2,k+1}^{n} + \\frac{\\Delta t}{2(\\Delta \\hat{z})^2}\\frac{\\mathcal{K}_{\\mathrm{M},j+1/2,k-1/2}}{J_{j+1/2}}v_{j+1/2,k-1}^{n} .
\\end{align*}
```

```julia 
turbulent_diffusion!(state::State, dt::AbstractFloat, variable::W)
```

Apply diffusion to the vertical momentum. 

The prognostic equation

```math
\\frac{\\partial w}{\\partial t} = \\frac{1}{J}\\frac{\\partial}{\\partial \\hat{z}}\\left(\\frac{K_\\mathrm{M}}{J}\\frac{\\partial w}{\\partial \\hat{z}}\\right)
```

is solved using the Crank-Nicolson scheme, where the system 

```math
a_{k+1/2} w_{k-1/2}^{n+1} + b_{k+1/2} w_{k+1/2}^{n+1} + c_{k+1/2} w_{k+3/2}^{n+1} = f_{k+1/2}
```

is solved using a Thomas tridiagonal solver, with ``\\mathcal{K}_\\mathrm{M} = \\frac{K_\\mathrm{M}}{J}`` and

```math 
\\begin{align*}
    a_{k+1/2} &= -\\frac{\\Delta t}{2(\\Delta \\hat{z})^2}\\frac{\\mathcal{K}_{\\mathrm{M},k}}{J_{k+1/2}^2} ,\\\\
    b_{k+1/2} &= 1 + \\frac{\\Delta t}{2(\\Delta \\hat{z})^2}\\frac{\\mathcal{K}_{\\mathrm{M},k+1}}{J_{k+1/2}^2} 
        + \\frac{\\Delta t}{2(\\Delta \\hat{z})^2}\\frac{\\mathcal{K}_{\\mathrm{M},k}}{J_{k+1/2}^2} ,\\\\
    c_{k+1/2} &= -\\frac{\\Delta t}{2(\\Delta \\hat{z})^2}\\frac{\\mathcal{K}_{\\mathrm{M},k+1}}{J_{k+1/2}^2} ,\\\\
    f_{k+1/2} &= \\left[1 - \\frac{\\Delta t}{2(\\Delta \\hat{z})^2}\\frac{\\mathcal{K}_{\\mathrm{M},k+1}}{J_{k+1/2}} 
        - \\frac{\\Delta t}{2(\\Delta \\hat{z})^2}\\frac{\\mathcal{K}_{\\mathrm{M},k}}{J_{k+1/2}}\\right]w_{k+1/2}^{n}\\\\
        &+\\frac{\\Delta t}{2(\\Delta \\hat{z})^2}\\frac{\\mathcal{K}_{\\mathrm{M},k+1}}{J_{k+1/2}}w_{k+1}^{n} + \\frac{\\Delta t}{2(\\Delta \\hat{z})^2}\\frac{\\mathcal{K}_{\\mathrm{M},k}}{J_{k+1/2}}w_{k-1}^{n} .
\\end{align*}
```

Using the relation 

```math
\\frac{\\partial \\hat{w}}{\\partial t} = G^{13}\\frac{\\partial u}{\\partial t} + G^{23}\\frac{\\partial v}{\\partial t} + \\frac{1}{J}\\frac{\\partial w}{\\partial t} ,
```

the transformed wind is calculated:

```math 
\\begin{align*}
\\hat{w}^{n+1}_{k+1/2} &= G_{k+1/2}^{13}u^{n+1}_{k+1/2} + G_{k+1/2}^{23}v^{n+1}_{k+1/2} + \\frac{1}{J_{k+1/2}}w^{n+1}_{k+1/2}.
\\end{align*}
```


```julia 
turbulent_diffusion!(state::State, dt::AbstractFloat, variable::Theta)
```

Apply diffusion to the mass-weighted potential temperature by dispatching to model-specific methods.

```julia 
turbulent_diffusion!(
    state::State,
    dt::AbstractFloat,
    variable::Theta,
    model::Union{PseudoIncompressible, Boussinesq},
)
```

Return for configurations in Boussinesq and pseudo-incompressible mode.

```julia 
turbulent_diffusion!(
    state::State,
    dt::AbstractFloat,
    variable::Theta,
    model::Compressible,
)
```

Apply diffusion to the potential temperature for configurations in Compressible mode.

The prognostic equation 

```math 
\\frac{\\partial \\theta}{\\partial t} = \\frac{1}{J}\\frac{\\partial}{\\partial \\hat{z}}\\left(\\frac{K_\\mathrm{H}}{J}\\frac{\\partial \\theta}{\\partial \\hat{z}}\\right)
```

is solved using the Crank-Nicolson scheme, where the system 

```math
a_k \\theta_{k-1}^{n+1} + b_k \\theta_k^{n+1} + c_k \\theta_{k+1}^{n+1} = f_k
```

is solved using a Thomas tridiagonal solver, with ``\\mathcal{K}_\\mathrm{H} = \\frac{K_\\mathrm{H}}{J}`` and 

```math 
\\begin{align*}
    a_k &= -\\frac{\\Delta t}{2(\\Delta \\hat{z})^2}\\frac{\\mathcal{K}_{\\mathrm{H},k-1/2}}{J_k}  , \\\\
    b_k &= 1 + \\frac{\\Delta t}{2(\\Delta \\hat{z})^2}\\frac{\\mathcal{K}_{\\mathrm{H},k+1/2}}{J_k} + \\frac{\\Delta t}{2(\\Delta \\hat{z})^2}\\frac{\\mathcal{K}_{\\mathrm{H},k-1/2}{J_k}}, \\\\
    c_k &= -\\frac{\\Delta t}{2(\\Delta \\hat{z})^2}\\frac{\\mathcal{K}_{\\mathrm{H},k+1/2}}{J_k}  , \\\\
    f_k &= \\left[1 - \\frac{\\Delta t}{2(\\Delta \\hat{z})^2}\\frac{\\mathcal{K}_{\\mathrm{H},k+1/2}}{J_k}  - \\frac{\\Delta t}{2(\\Delta \\hat{z})^2}\\frac{\\mathcal{K}_{\\mathrm{H},k-1/2}}{J_k}\\right] \\left(\\theta\\right)_k^{n} \\\\
    & + \\frac{\\Delta t}{2(\\Delta \\hat{z})^2}\\frac{\\mathcal{K}_{\\mathrm{H},k+1/2}}{J_k} \\left(\\theta\\right)_{k+1}^{n} + \\frac{\\Delta t}{2(\\Delta \\hat{z})^2}\\frac{\\mathcal{K}_{\\mathrm{H},k-1/2}}{J_k}  \\left(\\theta\\right)_{k-1}^{n}.
\\end{align*}
```

```julia 
turbulent_diffusion!(state::State, dt::AbstractFloat, variable::Chi)
```

Apply diffusion to tracers by dispatching to tracer-setup-specific configurations.

```julia 
turbulent_diffusion!(
    state::State,
    dt::AbstractFloat,
    variable::Chi,
    tracer_setup::NoTracer,
)
```

Return for configurations without tracer transport.

```julia 
turbulent_diffusion!(
    state::State,
    dt::AbstractFloat,
    variable::Chi,
    tracer_setup::TracerOn,
)
```

Apply diffusion to the tracers variables.

The prognostic equation 

```math 
\\frac{\\partial \\chi}{\\partial t} = \\frac{1}{J}\\frac{\\partial}{\\partial \\hat{z}}\\left(\\frac{K_\\mathrm{H}}{J}\\frac{\\partial \\chi}{\\partial \\hat{z}}\\right)
```

is solved using the Crank-Nicolson scheme, where the system 

```math
a_k \\chi_{k-1}^{n+1} + b_k \\chi_k^{n+1} + c_k \\chi_{k+1}^{n+1} = f_k
```

is solved using a Thomas tridiagonal solver, with ``\\mathcal{K}_\\mathrm{H} = \\frac{K_\\mathrm{H}}{J}`` and 

```math 
\\begin{align*}
    a_k &= -\\frac{\\Delta t}{2(\\Delta \\hat{z})^2}\\frac{\\mathcal{K}_{\\mathrm{H},k-1/2}}{J_k}  , \\\\
    b_k &= 1 + \\frac{\\Delta t}{2(\\Delta \\hat{z})^2}\\frac{\\mathcal{K}_{\\mathrm{H},k+1/2}}{J_k} + \\frac{\\Delta t}{2(\\Delta \\hat{z})^2}\\frac{\\mathcal{K}_{\\mathrm{H},k-1/2}}{J_k}, \\\\
    c_k &= -\\frac{\\Delta t}{2(\\Delta \\hat{z})^2}\\frac{\\mathcal{K}_{\\mathrm{H},k+1/2}}{J_k}  , \\\\
    f_k &= \\left[1 - \\frac{\\Delta t}{2(\\Delta \\hat{z})^2}\\frac{\\mathcal{K}_{\\mathrm{H},k+1/2}}{J_k}  - \\frac{\\Delta t}{2(\\Delta \\hat{z})^2}\\frac{\\mathcal{K}_{\\mathrm{H},k-1/2}}{J_k}\\right] \\chi_k^{n} \\\\
    &+ \\frac{\\Delta t}{2(\\Delta \\hat{z})^2}\\frac{\\mathcal{K}_{\\mathrm{H},k+1/2}}{J_k} \\chi_{k+1}^{n} + \\frac{\\Delta t}{2(\\Delta \\hat{z})^2}\\frac{\\mathcal{K}_{\\mathrm{H},k-1/2}}{J_k}  \\chi_{k-1}^{n}.
\\end{align*}
```

# Arguments

  - `state`: Model state. 

  - `dt`: Time step.

  - `turbulence_scheme`: General turbulence parameterization configuration.

  - `variable`: Variable (equation) of choice.

  - `model`: Dynamic equations.

  - `tracer_setup`: General tracer-transport configuration.

# See also 

  - [`PinCFlow.Update.check_tke!`](@ref)

  - [`PinCFlow.Update.reset_thomas!`](@ref)

  - [`PinCFlow.Update.thomas_algorithm!`](@ref)

  - [`PinCFlow.Update.turbulence_diffusion_coefficient`](@ref)
"""
function turbulent_diffusion! end

function turbulent_diffusion!(state::State, dt::AbstractFloat)
    (; turbulence_scheme) = state.namelists.turbulence

    @dispatch_turbulence_scheme turbulent_diffusion!(
        state,
        dt,
        Val(turbulence_scheme),
    )
    return
end

function turbulent_diffusion!(
    state::State,
    dt::AbstractFloat,
    turbulence_scheme::Val{:NoTurbulence},
)
    return
end

function turbulent_diffusion!(
    state::State,
    dt::AbstractFloat,
    turbulence_scheme::Val{:TKEScheme},
)
    (; momentum_coupling, entropy_coupling, tracer_coupling) =
        state.namelists.turbulence

    check_tke!(state)

    if momentum_coupling
        turbulent_diffusion!(state, dt, U())
        turbulent_diffusion!(state, dt, V())
        turbulent_diffusion!(state, dt, W())
    end
    if entropy_coupling
        turbulent_diffusion!(state, dt, Theta())
    end
    if tracer_coupling
        turbulent_diffusion!(state, dt, Chi())
    end
    return
end

function turbulent_diffusion!(state::State, dt::AbstractFloat, variable::U)
    (; u) = state.variables.predictands
    (; i0, i1, j0, j1, k0, k1) = state.domain
    (; jac, dz) = state.grid
    (; ath, bth, cth, fth) = state.variables.auxiliaries

    dtdz2 = dt / (2.0 * dz^2.0)

    reset_thomas!(state)

    @ivy for k in k0:k1, j in j0:j1, i in i0:i1
        kmu =
            0.5 * (
                (
                    jac[i, j, k] * (
                        turbulence_diffusion_coefficient(
                            state,
                            i,
                            j,
                            k + 1,
                            KM(),
                        ) / jac[i, j, k + 1]
                    ) +
                    jac[i, j, k + 1] * (
                        turbulence_diffusion_coefficient(state, i, j, k, KM()) /
                        jac[i, j, k]
                    )
                ) / (jac[i, j, k] + jac[i, j, k + 1]) +
                (
                    jac[i + 1, j, k] * (
                        turbulence_diffusion_coefficient(
                            state,
                            i + 1,
                            j,
                            k + 1,
                            KM(),
                        ) / jac[i + 1, j, k + 1]
                    ) +
                    jac[i + 1, j, k + 1] * (
                        turbulence_diffusion_coefficient(
                            state,
                            i + 1,
                            j,
                            k,
                            KM(),
                        ) / jac[i + 1, j, k]
                    )
                ) / (jac[i + 1, j, k] + jac[i + 1, j, k + 1])
            )

        kmd =
            0.5 * (
                (
                    jac[i, j, k] * (
                        turbulence_diffusion_coefficient(
                            state,
                            i,
                            j,
                            k + 1,
                            KM(),
                        ) / jac[i, j, k - 1]
                    ) +
                    jac[i, j, k - 1] * (
                        turbulence_diffusion_coefficient(state, i, j, k, KM()) /
                        jac[i, j, k]
                    )
                ) / (jac[i, j, k] + jac[i, j, k - 1]) +
                (
                    jac[i + 1, j, k] * (
                        turbulence_diffusion_coefficient(
                            state,
                            i + 1,
                            j,
                            k - 1,
                            KM(),
                        ) / jac[i + 1, j, k - 1]
                    ) +
                    jac[i + 1, j, k - 1] * (
                        turbulence_diffusion_coefficient(
                            state,
                            i + 1,
                            j,
                            k,
                            KM(),
                        ) / jac[i + 1, j, k]
                    )
                ) / (jac[i + 1, j, k] + jac[i + 1, j, k - 1])
            )

        jacc = 0.5 * (jac[i, j, k] + jac[i + 1, j, k])

        ith = i - i0 + 1
        jth = j - j0 + 1
        kth = k - k0 + 1

        ath[ith, jth, kth] = -dtdz2 * kmd / jacc
        bth[ith, jth, kth] = 1 + dtdz2 * kmu / jacc + dtdz2 * kmd / jacc
        cth[ith, jth, kth] = -dtdz2 * kmu / jacc
        fth[ith, jth, kth] =
            dtdz2 * kmu / jacc * u[i, j, k + 1] +
            (1 - dtdz2 * kmu / jacc - dtdz2 * kmd / jacc) * u[i, j, k] +
            dtdz2 * kmd / jacc * u[i, j, k - 1]
    end

    thomas_algorithm!(state)

    @ivy u[i0:i1, j0:j1, k0:k1] .= fth
    return
end

function turbulent_diffusion!(state::State, dt::AbstractFloat, variable::V)
    (; v) = state.variables.predictands
    (; i0, i1, j0, j1, k0, k1) = state.domain
    (; jac, dz) = state.grid
    (; ath, bth, cth, fth) = state.variables.auxiliaries

    dtdz2 = dt / (2.0 * dz^2.0)

    reset_thomas!(state)

    @ivy for k in k0:k1, j in j0:j1, i in i0:i1
        kmu =
            0.5 * (
                (
                    jac[i, j, k] * (
                        turbulence_diffusion_coefficient(
                            state,
                            i,
                            j,
                            k + 1,
                            KM(),
                        ) / jac[i, j, k + 1]
                    ) +
                    jac[i, j, k + 1] * (
                        turbulence_diffusion_coefficient(state, i, j, k, KM()) /
                        jac[i, j, k]
                    )
                ) / (jac[i, j, k] + jac[i, j, k + 1]) +
                (
                    jac[i, j + 1, k] * (
                        turbulence_diffusion_coefficient(
                            state,
                            i,
                            j + 1,
                            k + 1,
                            KM(),
                        ) / jac[i, j + 1, k + 1]
                    ) +
                    jac[i, j + 1, k + 1] * (
                        turbulence_diffusion_coefficient(
                            state,
                            i,
                            j + 1,
                            k,
                            KM(),
                        ) / jac[i, j + 1, k]
                    )
                ) / (jac[i, j + 1, k] + jac[i, j + 1, k + 1])
            )

        kmd =
            0.5 * (
                (
                    jac[i, j, k] * (
                        turbulence_diffusion_coefficient(
                            state,
                            i,
                            j,
                            k - 1,
                            KM(),
                        ) / jac[i, j, k - 1]
                    ) +
                    jac[i, j, k - 1] * (
                        turbulence_diffusion_coefficient(state, i, j, k, KM()) /
                        jac[i, j, k]
                    )
                ) / (jac[i, j, k] + jac[i, j, k - 1]) +
                (
                    jac[i, j + 1, k] * (
                        turbulence_diffusion_coefficient(
                            state,
                            i,
                            j + 1,
                            k - 1,
                            KM(),
                        ) / jac[i, j + 1, k - 1]
                    ) +
                    jac[i, j + 1, k - 1] * (
                        turbulence_diffusion_coefficient(
                            state,
                            i,
                            j + 1,
                            k,
                            KM(),
                        ) / jac[i, j + 1, k]
                    )
                ) / (jac[i, j + 1, k] + jac[i, j + 1, k - 1])
            )

        jacc = 0.5 * (jac[i, j, k] + jac[i, j + 1, k])

        ith = i - i0 + 1
        jth = j - j0 + 1
        kth = k - k0 + 1

        ath[ith, jth, kth] = -dtdz2 * kmd / jacc
        bth[ith, jth, kth] = 1 + dtdz2 * kmu / jacc + dtdz2 * kmd / jacc
        cth[ith, jth, kth] = -dtdz2 * kmu / jacc
        fth[ith, jth, kth] =
            dtdz2 * kmu / jacc * v[i, j, k + 1] +
            (1 - dtdz2 * kmu / jacc - dtdz2 * kmd / jacc) * v[i, j, k] +
            dtdz2 * kmd / jacc * v[i, j, k - 1]
    end

    thomas_algorithm!(state)

    @ivy v[i0:i1, j0:j1, k0:k1] .= fth
    return
end

function turbulent_diffusion!(state::State, dt::AbstractFloat, variable::W)
    (; u, v, w) = state.variables.predictands
    (; i0, i1, j0, j1, k0, k1) = state.domain
    (; jac, met, dz) = state.grid
    (; ath, bth, cth, fth) = state.variables.auxiliaries

    dtdz2 = dt / (2.0 * dz^2.0)

    reset_thomas!(state)

    @ivy for k in k0:k1, j in j0:j1, i in i0:i1
        kmu =
            turbulence_diffusion_coefficient(state, i, j, k + 1, KM()) /
            jac[i, j, k + 1]
        kmd =
            turbulence_diffusion_coefficient(state, i, j, k - 1, KM()) /
            jac[i, j, k - 1]

        wu = compute_vertical_wind(i, j, k + 1, state)
        wc = compute_vertical_wind(i, j, k, state)
        wd = compute_vertical_wind(i, j, k - 1, state)

        jacc =
            2.0 * jac[i, j, k] * jac[i, j, k + 1] /
            (jac[i, j, k] + jac[i, j, k + 1])

        ith = i - i0 + 1
        jth = j - j0 + 1
        kth = k - k0 + 1

        ath[ith, jth, kth] = -dtdz2 / jacc * kmd
        bth[ith, jth, kth] = 1 + dtdz2 / jacc * kmu + dtdz2 / jacc^2 * kmd
        cth[ith, jth, kth] = -dtdz2 / jacc * kmu
        fth[ith, jth, kth] =
            dtdz2 / jacc * kmu * wu +
            (1 - dtdz2 / jacc * kmu - dtdz2 / jacc^2 * kmd) * wc +
            dtdz2 / jacc * kmd * wd
    end

    thomas_algorithm!(state)

    @ivy for k in k0:k1, j in j0:j1, i in i0:i1
        ud = 0.5 * (u[i - 1, j, k] + u[i, j, k])
        uu = 0.5 * (u[i - 1, j, k + 1] + u[i, j, k + 1])
        uc13 =
            (
                jac[i, j, k] * met[i, j, k + 1, 1, 3] * uu +
                jac[i, j, k + 1] * met[i, j, k, 1, 3] * ud
            ) / (jac[i, j, k] + jac[i, j, k + 1])

        vd = 0.5 * (v[i, j - 1, k] + v[i, j, k])
        vu = 0.5 * (v[i, j - 1, k + 1] + v[i, j, k + 1])
        vc23 =
            (
                jac[i, j, k] * met[i, j, k + 1, 2, 3] * vu +
                jac[i, j, k + 1] * met[i, j, k, 2, 3] * vd
            ) / (jac[i, j, k] + jac[i, j, k + 1])

        jacc =
            2.0 * jac[i, j, k] * jac[i, j, k + 1] /
            (jac[i, j, k] + jac[i, j, k + 1])

        ith = i - i0 + 1
        jth = j - j0 + 1
        kth = k - k0 + 1

        w[i, j, k] = uc13 + vc23 + fth[ith, jth, kth] / jacc
    end

    return
end

function turbulent_diffusion!(state::State, dt::AbstractFloat, variable::Theta)
    (; model) = state.namelists.atmosphere

    @dispatch_model turbulent_diffusion!(state, dt, variable, Val(model))
    return
end

function turbulent_diffusion!(
    state::State,
    dt::AbstractFloat,
    variable::Theta,
    model::Union{Val{:PseudoIncompressible}, Val{:Boussinesq}},
)
    return
end

function turbulent_diffusion!(
    state::State,
    dt::AbstractFloat,
    variable::Theta,
    model::Val{:Compressible},
)
    (; p, rho) = state.variables.predictands
    (; i0, i1, j0, j1, k0, k1) = state.domain
    (; jac, dz) = state.grid
    (; ath, bth, cth, fth) = state.variables.auxiliaries
    (; rhobar) = state.atmosphere

    dtdz2 = dt / (2.0 * dz^2.0)

    reset_thomas!(state)

    @ivy for k in k0:k1, j in j0:j1, i in i0:i1
        khd =
            (
                jac[i, j, k - 1] * (
                    turbulence_diffusion_coefficient(state, i, j, k, KH()) /
                    jac[i, j, k]
                ) +
                jac[i, j, k] * (
                    turbulence_diffusion_coefficient(state, i, j, k - 1, KH()) /
                    jac[i, j, k - 1]
                )
            ) / (jac[i, j, k - 1] + jac[i, j, k])
        khu =
            (
                jac[i, j, k + 1] * (
                    turbulence_diffusion_coefficient(state, i, j, k, KH()) /
                    jac[i, j, k]
                ) +
                jac[i, j, k] * (
                    turbulence_diffusion_coefficient(state, i, j, k + 1, KH()) /
                    jac[i, j, k + 1]
                )
            ) / (jac[i, j, k + 1] + jac[i, j, k])

        ith = i - i0 + 1
        jth = j - j0 + 1
        kth = k - k0 + 1

        ath[ith, jth, kth] = -dtdz2 / jac[i, j, k] * khd
        bth[ith, jth, kth] =
            1 + dtdz2 / jac[i, j, k] * khu + dtdz2 / jac[i, j, k] * khd
        cth[ith, jth, kth] = -dtdz2 / jac[i, j, k] * khu

        fth[ith, jth, kth] =
            (1 - dtdz2 / jac[i, j, k] * khu - dtdz2 / jac[i, j, k] * khd) *
            p[i, j, k] / (rho[i, j, k] + rhobar[i, j, k]) +
            dtdz2 / jac[i, j, k] * khu * p[i, j, k + 1] /
            (rho[i, j, k + 1] + rhobar[i, j, k + 1]) +
            dtdz2 / jac[i, j, k] * khd * p[i, j, k - 1] /
            (rho[i, j, k - 1] + rhobar[i, j, k - 1])
    end

    thomas_algorithm!(state)

    @ivy p[i0:i1, j0:j1, k0:k1] .=
        fth .* (rho[i0:i1, j0:j1, k0:k1] .+ rhobar[i0:i1, j0:j1, k0:k1])
    return
end

function turbulent_diffusion!(state::State, dt::AbstractFloat, variable::Chi)
    (; tracer_setup) = state.namelists.tracer

    @dispatch_tracer_setup turbulent_diffusion!(
        state,
        dt,
        variable,
        Val(tracer_setup),
    )
    return
end

function turbulent_diffusion!(
    state::State,
    dt::AbstractFloat,
    variable::Chi,
    tracer_setup::Val{:NoTracer},
)
    return
end

function turbulent_diffusion!(
    state::State,
    dt::AbstractFloat,
    variable::Chi,
    tracer_setup::Val{:TracerOn},
)
    (; tracerpredictands) = state.tracer
    (; rho) = state.variables.predictands
    (; rhobar) = state.atmosphere
    (; i0, i1, j0, j1, k0, k1) = state.domain
    (; jac, dz) = state.grid
    (; ath, bth, cth, fth) = state.variables.auxiliaries

    dtdz2 = dt / (2.0 * dz^2.0)

    reset_thomas!(state)

    @ivy for field in 1:fieldcount(TracerPredictands)
        chi = getfield(tracerpredictands, field)
        for k in k0:k1, j in j0:j1, i in i0:i1
            khd =
                (
                    jac[i, j, k - 1] * (
                        turbulence_diffusion_coefficient(state, i, j, k, KH()) /
                        jac[i, j, k]
                    ) +
                    jac[i, j, k] * (
                        turbulence_diffusion_coefficient(
                            state,
                            i,
                            j,
                            k - 1,
                            KH(),
                        ) / jac[i, j, k - 1]
                    )
                ) / (jac[i, j, k - 1] + jac[i, j, k])
            khu =
                (
                    jac[i, j, k + 1] * (
                        turbulence_diffusion_coefficient(state, i, j, k, KH()) /
                        jac[i, j, k]
                    ) +
                    jac[i, j, k] * (
                        turbulence_diffusion_coefficient(
                            state,
                            i,
                            j,
                            k + 1,
                            KH(),
                        ) / jac[i, j, k + 1]
                    )
                ) / (jac[i, j, k + 1] + jac[i, j, k])

            ith = i - i0 + 1
            jth = j - j0 + 1
            kth = k - k0 + 1

            ath[ith, jth, kth] = -dtdz2 / jac[i, j, k] * khd
            bth[ith, jth, kth] =
                1 + dtdz2 / jac[i, j, k] * khu + dtdz2 / jac[i, j, k] * khd
            cth[ith, jth, kth] = -dtdz2 / jac[i, j, k] * khu

            fth[ith, jth, kth] =
                (1 - dtdz2 / jac[i, j, k] * khu - dtdz2 / jac[i, j, k] * khd) *
                chi[i, j, k] / (rho[i, j, k] + rhobar[i, j, k]) +
                dtdz2 / jac[i, j, k] * khu * chi[i, j, k + 1] /
                (rho[i, j, k + 1] + rhobar[i, j, k + 1]) +
                dtdz2 / jac[i, j, k] * khd * chi[i, j, k - 1] /
                (rho[i, j, k - 1] + rhobar[i, j, k - 1])
        end

        thomas_algorithm!(state)

        chi[i0:i1, j0:j1, k0:k1] .=
            fth .* (rho[i0:i1, j0:j1, k0:k1] .+ rhobar[i0:i1, j0:j1, k0:k1])
    end
    return
end
