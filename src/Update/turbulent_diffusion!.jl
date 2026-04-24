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
\\frac{\\partial u}{\\partial t} = \\frac{1}{J}\\frac{\\partial}{\\partial \\hat{z}}\\left(J K_\\mathrm{M}G^{33}\\frac{\\partial u}{\\partial \\hat{z}}\\right)
```

is solved using the Crank-Nicolson scheme, where the system 

```math
a_{i+1/2,k} u_{i+1/2,k-1}^{n+1} + b_{i+1/2,k} u_{i+1/2,k}^{n+1} + c_{i+1/2,k} u_{i+1/2,k+1}^{n+1} = f_{i+1/2,k}
```

is solved using a Thomas tridiagonal solver, with ``\\mathcal{K}^{33} = J K_\\mathrm{M}G^{33}`` and

```math 
\\begin{align*}
    a_{i+1/2,k} =& -\\frac{\\Delta t}{2(\\Delta \\hat{z})^2}\\frac{\\mathcal{K}^{33}_{\\mathrm{M},i+1/2,k-1/2}}{J_{i+1/2}} ,\\\\
    b_{i+1/2,k} =& \\left(1 + \\frac{\\Delta t}{2(\\Delta \\hat{z})^2}\\frac{\\mathcal{K}^{33}_{\\mathrm{M},i+1/2,k+1/2}}{J_{i+1/2}} 
        + \\frac{\\Delta t}{2(\\Delta \\hat{z})^2}\\frac{\\mathcal{K}^{33}_{\\mathrm{M},i+1/2,k-1/2}}{J_{i+1/2}}\\right) ,\\\\
    c_{i+1/2,k} =& -\\frac{\\Delta t}{2(\\Delta \\hat{z})^2}\\frac{\\mathcal{K}^{33}_{\\mathrm{M},i+1/2,k+1/2}}{J_{i+1/2}} ,\\\\
    f_{i+1/2,k} =& \\left(1 - \\frac{\\Delta t}{2(\\Delta \\hat{z})^2}\\frac{\\mathcal{K}^{33}_{\\mathrm{M},i+1/2,k+1/2}}{J_{i+1/2}} 
        - \\frac{\\Delta t}{2(\\Delta \\hat{z})^2}\\frac{\\mathcal{K}^{33}_{\\mathrm{M},i+1/2,k-1/2}}{J_{i+1/2}}\\right)v_{i+1/2}^{n} \\\\
        &+\\frac{\\Delta t}{2(\\Delta \\hat{z})^2}\\frac{\\mathcal{K}^{33}_{\\mathrm{M},i+1/2,k+1/2}}{J_{i+1/2}}v_{i+1/2,k+1}^{n}
        +\\frac{\\Delta t}{2(\\Delta \\hat{z})^2}\\frac{\\mathcal{K}^{33}_{\\mathrm{M},i+1/2,k-1/2}}{J_{i+1/2}}v_{i+1/2,k-1}^{n} .
\\end{align*}
```

```julia 
turbulent_diffusion!(state::State, dt::AbstractFloat, variable::V)
```

Apply diffusion to the meridional momentum. 

The prognostic equation

```math
\\frac{\\partial v}{\\partial t} = \\frac{1}{J}\\frac{\\partial}{\\partial \\hat{z}}\\left(J K_\\mathrm{M}G^{33}\\frac{\\partial v}{\\partial \\hat{z}}\\right)
```

is solved using the Crank-Nicolson scheme, where the system 

```math
a_{j+1/2,k} v_{j+1/2,k-1}^{n+1} + b_{j+1/2,k} v_{j+1/2,k}^{n+1} + c_{j+1/2,k} v_{j+1/2,k+1}^{n+1} = f_{j+1/2,k}
```

is solved using a Thomas tridiagonal solver, with ``\\mathcal{K}^{33} = J K_\\mathrm{M}G^{33}`` and

```math 
\\begin{align*}
    a_{j+1/2,k} =& -\\frac{\\Delta t}{2(\\Delta \\hat{z})^2}\\frac{\\mathcal{K}^{33}_{\\mathrm{M},j+1/2,k-1/2}}{J_{j+1/2}} ,\\\\
    b_{j+1/2,k} =& \\left(1 + \\frac{\\Delta t}{2(\\Delta \\hat{z})^2}\\frac{\\mathcal{K}^{33}_{\\mathrm{M},j+1/2,k+1/2}}{J_{j+1/2}} 
        + \\frac{\\Delta t}{2(\\Delta \\hat{z})^2}\\frac{\\mathcal{K}^{33}_{\\mathrm{M},j+1/2,k-1/2}}{J_{j+1/2}}\\right) ,\\\\
    c_{j+1/2,k} =& -\\frac{\\Delta t}{2(\\Delta \\hat{z})^2}\\frac{\\mathcal{K}^{33}_{\\mathrm{M},j+1/2,k+1/2}}{J_{j+1/2}} ,\\\\
    f_{j+1/2,k} =& \\left(1 - \\frac{\\Delta t}{2(\\Delta \\hat{z})^2}\\frac{\\mathcal{K}^{33}_{\\mathrm{M},j+1/2,k+1/2}}{J_{j+1/2}} 
        - \\frac{\\Delta t}{2(\\Delta \\hat{z})^2}\\frac{\\mathcal{K}^{33}_{\\mathrm{M},j+1/2,k-1/2}}{J_{j+1/2}}\\right)v_{j+1/2}^{n}\\\\
        &+\\frac{\\Delta t}{2(\\Delta \\hat{z})^2}\\frac{\\mathcal{K}^{33}_{\\mathrm{M},j+1/2,k+1/2}}{J_{j+1/2}}v_{j+1/2,k+1}^{n}
        &+\\frac{\\Delta t}{2(\\Delta \\hat{z})^2}\\frac{\\mathcal{K}^{33}_{\\mathrm{M},j+1/2,k-1/2}}{J_{j+1/2}}v_{j+1/2,k-1}^{n} .
\\end{align*}
```

```julia 
turbulent_diffusion!(state::State, dt::AbstractFloat, variable::W)
```

Apply diffusion to the vertical momentum. 

The prognostic equation

```math
\\frac{\\partial \\hat{w}}{\\partial t} = G^{13}\\frac{\\partial u}{\\partial t} + G^{23}\\frac{\\partial v}{\\partial t} + \\frac{1}{J}\\frac{\\partial}{\\partial \\hat{z}}\\left(K_\\mathrm{M}\\frac{\\partial w}{\\partial \\hat{z}}\\right)
```

is solved using the Crank-Nicolson scheme, where the system 

```math
a_{k+1/2} \\hat{w}_{k-1/2}^{n+1} + b_{k+1/2} \\hat{w}_{k+1/2}^{n+1} + c_{k+1/2} \\hat{w}_{k+3/2}^{n+1} = f_{k+1/2}
```

is solved using a Thomas tridiagonal solver, with ``\\mathcal{K}^{33} = J K_\\mathrm{M}G^{33}`` and

```math 
\\begin{align*}
    a_{k+1/2} =& -\\frac{\\Delta t}{2(\\Delta \\hat{z})^2}\\frac{\\mathcal{K}^{33}_{\\mathrm{M},k}}{J_{k+1/2}} ,\\\\
    b_{k+1/2} =& \\left(1 + \\frac{\\Delta t}{2(\\Delta \\hat{z})^2}\\frac{\\mathcal{K}^{33}_{\\mathrm{M},k+1}}{J_{k+1/2}} 
        + \\frac{\\Delta t}{2(\\Delta \\hat{z})^2}\\frac{\\mathcal{K}^{33}_{\\mathrm{M},k}}{J_{k+1/2}}\\right) ,\\\\
    c_{k+1/2} =& -\\frac{\\Delta t}{2(\\Delta \\hat{z})^2}\\frac{\\mathcal{K}^{33}_{\\mathrm{M},k+1}}{J_{k+1/2}} ,\\\\
    f_{k+1/2} =& G_{k+1/2}^{13}\\left(\\frac{u^{n+1}-u^{n}}{\\Delta t}\\right)_{k+1/2} + G_{k+1/2}^{23}\\left(\\frac{v^{n+1}-v^{n}}{\\Delta t}\\right)_{k+1/2} + \\left(1 - \\frac{\\Delta t}{2(\\Delta \\hat{z})^2}\\frac{\\mathcal{K}^{33}_{\\mathrm{M},k+1}}{J_{k+1/2}} 
        - \\frac{\\Delta t}{2(\\Delta \\hat{z})^2}\\frac{\\mathcal{K}^{33}_{\\mathrm{M},k}}{J_{k+1/2}}\\right)\\hat{w}_{k+1/2}^{n}\\\\
        &+\\frac{\\Delta t}{2(\\Delta \\hat{z})^2}\\frac{\\mathcal{K}^{33}_{\\mathrm{M},k+1}}{J_{k+1/2}}\\hat{w}_{k+1}^{n}
        &+\\frac{\\Delta t}{2(\\Delta \\hat{z})^2}\\frac{\\mathcal{K}^{33}_{\\mathrm{M},k}}{J_{k+1/2}}\\hat{w}_{k-1}^{n} .
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

Apply diffusion to the mass-weighted potential temperature for configurations in Compressible mode.

The prognostic equation 

```math 
\\frac{\\partial P}{\\partial t} = \\frac{1}{J}\\frac{\\partial}{\\partial \\hat{z}}\\left[J K_H G^{33}\\frac{\\partial P}{\\partial \\hat{z}}\\right]
```

is solved using the Crank-Nicolson scheme, where the system 

```math
a_k P_{k-1}^{n+1} + b_k P_k^{n+1} + c_k P_{k+1}^{n+1} = f_k
```

is solved using a Thomas tridiagonal solver, with ``\\mathcal{K}^{33} = J K_\\mathrm{H}G^{33}`` and 

```math 
\\begin{align*}
    a_k = & -\\frac{\\Delta t}{2(\\Delta \\hat{z})^2}\\frac{\\mathcal{K}_{\\mathrm{H},k-1/2}^{33}}{J_k}  , \\\\
    b_k = & 1 + \\frac{\\Delta t}{2(\\Delta \\hat{z})^2}\\frac{\\mathcal{K}_{\\mathrm{H},k+1/2}^{33}}{J_k} + \\frac{\\Delta t}{2(\\Delta \\hat{z})^2}\\frac{\\mathcal{K}_{\\mathrm{H},k-1/2}^{33}{J_k}}, \\\\
    c_k = & -\\frac{\\Delta t}{2(\\Delta \\hat{z})^2}\\frac{\\mathcal{K}_{\\mathrm{H},k+1/2}^{33}}{J_k}  , \\\\
    f_k = & \\left( 1 - \\frac{\\Delta t}{2(\\Delta \\hat{z})^2}\\frac{\\mathcal{K}_{\\mathrm{H},k+1/2}^{33}}{J_k}  - \\frac{\\Delta t}{2(\\Delta \\hat{z})^2}\\frac{\\mathcal{K}_{\\mathrm{H},k-1/2}^{33}}{J_k}\\right) \\left(\\rho\\chi\\right)_k^{n} \\\\
    & + \\frac{\\Delta t}{2(\\Delta \\hat{z})^2}\\frac{\\mathcal{K}_{\\mathrm{H},k+1/2}^{33}}{J_k} \\left(\\rho\\chi\\right)_{k+1}^{n} + \\frac{\\Delta t}{2(\\Delta \\hat{z})^2}\\frac{\\mathcal{K}_{\\mathrm{H},k-1/2}^{33}}{J_k}  \\left(\\rho\\chi\\right)_{k-1}^{n}.
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
\\frac{\\partial \\rho\\chi}{\\partial t} = \\frac{1}{J}\\frac{\\partial}{\\partial \\hat{z}}\\left[J K_H G^{33}\\frac{\\partial \\rho\\chi}{\\partial \\hat{z}}\\right]
```

is solved using the Crank-Nicolson scheme, where the system 

```math
a_k \\left(\\rho\\chi\\right)_{k-1}^{n+1} + b_k \\left(\\rho\\chi\\right)_k^{n+1} + c_k \\left(\\rho\\chi\\right)_{k+1}^{n+1} = f_k
```

is solved using a Thomas tridiagonal solver, with ``\\mathcal{K}^{33} = J K_\\mathrm{H}G^{33}`` and 

```math 
\\begin{align*}
    a_k = & -\\frac{\\Delta t}{2(\\Delta \\hat{z})^2}\\frac{\\mathcal{K}_{\\mathrm{H},k-1/2}^{33}}{J_k}  , \\\\
    b_k = & 1 + \\frac{\\Delta t}{2(\\Delta \\hat{z})^2}\\frac{\\mathcal{K}_{\\mathrm{H},k+1/2}^{33}}{J_k} + \\frac{\\Delta t}{2(\\Delta \\hat{z})^2}\\frac{\\mathcal{K}_{\\mathrm{H},k-1/2}^{33}{J_k}}, \\\\
    c_k = & -\\frac{\\Delta t}{2(\\Delta \\hat{z})^2}\\frac{\\mathcal{K}_{\\mathrm{H},k+1/2}^{33}}{J_k}  , \\\\
    f_k = & \\left( 1 - \\frac{\\Delta t}{2(\\Delta \\hat{z})^2}\\frac{\\mathcal{K}_{\\mathrm{H},k+1/2}^{33}}{J_k}  - \\frac{\\Delta t}{2(\\Delta \\hat{z})^2}\\frac{\\mathcal{K}_{\\mathrm{H},k-1/2}^{33}}{J_k}\\right) \\left(\\rho\\chi\\right)_k^{n} \\\\
    & + \\frac{\\Delta t}{2(\\Delta \\hat{z})^2}\\frac{\\mathcal{K}_{\\mathrm{H},k+1/2}^{33}}{J_k} \\left(\\rho\\chi\\right)_{k+1}^{n} + \\frac{\\Delta t}{2(\\Delta \\hat{z})^2}\\frac{\\mathcal{K}_{\\mathrm{H},k-1/2}^{33}}{J_k}  \\left(\\rho\\chi\\right)_{k-1}^{n}.
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
    (; uold, vold) = state.variables.backups
    (; u, v) = state.variables.predictands

    check_tke!(state)

    if momentum_coupling
        uold .= copy(u)
        vold .= copy(v)
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
    (; nbx, nby, nbz) = state.namelists.domain
    (; u) = state.variables.predictands
    (; i0, i1, j0, j1, k0, k1) = state.domain
    (; jac, met, dz) = state.grid
    (; ath, bth, cth, fth) = state.variables.auxiliaries

    dtdz2 = dt / (2.0 * dz^2.0)

    reset_thomas!(state)

    @ivy for k in k0:k1, j in j0:j1, i in i0:i1
        k33u =
            0.5 * (
                (
                    jac[i, j, k] * (
                        jac[i, j, k + 1] *
                        met[i, j, k + 1, 3, 3] *
                        turbulence_diffusion_coefficient(
                            state,
                            i,
                            j,
                            k + 1,
                            KM(),
                        )
                    ) +
                    jac[i, j, k + 1] * (
                        jac[i, j, k] *
                        met[i, j, k, 3, 3] *
                        turbulence_diffusion_coefficient(state, i, j, k, KM())
                    )
                ) / (jac[i, j, k] + jac[i, j, k + 1]) +
                (
                    jac[i + 1, j, k] * (
                        jac[i + 1, j, k + 1] *
                        met[i + 1, j, k + 1, 3, 3] *
                        turbulence_diffusion_coefficient(
                            state,
                            i + 1,
                            j,
                            k + 1,
                            KM(),
                        )
                    ) +
                    jac[i + 1, j, k + 1] * (
                        jac[i + 1, j, k] *
                        met[i + 1, j, k, 3, 3] *
                        turbulence_diffusion_coefficient(
                            state,
                            i + 1,
                            j,
                            k,
                            KM(),
                        )
                    )
                ) / (jac[i + 1, j, k] + jac[i + 1, j, k + 1])
            )

        k33d =
            0.5 * (
                (
                    jac[i, j, k] * (
                        jac[i, j, k - 1] *
                        met[i, j, k - 1, 3, 3] *
                        turbulence_diffusion_coefficient(
                            state,
                            i,
                            j,
                            k + 1,
                            KM(),
                        )
                    ) +
                    jac[i, j, k - 1] * (
                        jac[i, j, k] *
                        met[i, j, k, 3, 3] *
                        turbulence_diffusion_coefficient(state, i, j, k, KM())
                    )
                ) / (jac[i, j, k] + jac[i, j, k - 1]) +
                (
                    jac[i + 1, j, k] * (
                        jac[i + 1, j, k - 1] *
                        met[i + 1, j, k - 1, 3, 3] *
                        turbulence_diffusion_coefficient(
                            state,
                            i + 1,
                            j,
                            k - 1,
                            KM(),
                        )
                    ) +
                    jac[i + 1, j, k - 1] * (
                        jac[i + 1, j, k] *
                        met[i + 1, j, k, 3, 3] *
                        turbulence_diffusion_coefficient(
                            state,
                            i + 1,
                            j,
                            k,
                            KM(),
                        )
                    )
                ) / (jac[i + 1, j, k] + jac[i + 1, j, k - 1])
            )

        jacc = 0.5 * (jac[i, j, k] + jac[i + 1, j, k])

        ith = i - nbx
        jth = j - nby
        kth = k - nbz

        ath[ith, jth, kth] = -dtdz2 * k33d / jacc
        bth[ith, jth, kth] = 1 + dtdz2 * k33u / jacc + dtdz2 * k33d / jacc
        cth[ith, jth, kth] = -dtdz2 * k33u / jacc
        fth[ith, jth, kth] =
            dtdz2 * k33u / jacc * u[i, j, k + 1] +
            (1 - dtdz2 * k33u / jacc - dtdz2 * k33d / jacc) * u[i, j, k] +
            dtdz2 * k33d / jacc * u[i, j, k - 1]
    end

    thomas_algorithm!(state)

    u[i0:i1, j0:j1, k0:k1] .= fth
    return
end

function turbulent_diffusion!(state::State, dt::AbstractFloat, variable::V)
    (; nbx, nby, nbz) = state.namelists.domain
    (; v) = state.variables.predictands
    (; i0, i1, j0, j1, k0, k1) = state.domain
    (; jac, met, dz) = state.grid
    (; ath, bth, cth, fth) = state.variables.auxiliaries

    dtdz2 = dt / (2.0 * dz^2.0)

    reset_thomas!(state)

    @ivy for k in k0:k1, j in j0:j1, i in i0:i1
        k33u =
            0.5 * (
                (
                    jac[i, j, k] * (
                        jac[i, j, k + 1] *
                        met[i, j, k + 1, 3, 3] *
                        turbulence_diffusion_coefficient(
                            state,
                            i,
                            j,
                            k + 1,
                            KM(),
                        )
                    ) +
                    jac[i, j, k + 1] * (
                        jac[i, j, k] *
                        met[i, j, k, 3, 3] *
                        turbulence_diffusion_coefficient(state, i, j, k, KM())
                    )
                ) / (jac[i, j, k] + jac[i, j, k + 1]) +
                (
                    jac[i, j + 1, k] * (
                        jac[i, j + 1, k + 1] *
                        met[i, j + 1, k + 1, 3, 3] *
                        turbulence_diffusion_coefficient(
                            state,
                            i,
                            j + 1,
                            k + 1,
                            KM(),
                        )
                    ) +
                    jac[i, j + 1, k + 1] * (
                        jac[i, j + 1, k] *
                        met[i, j + 1, k, 3, 3] *
                        turbulence_diffusion_coefficient(
                            state,
                            i,
                            j + 1,
                            k,
                            KM(),
                        )
                    )
                ) / (jac[i, j + 1, k] + jac[i, j + 1, k + 1])
            )

        k33d =
            0.5 * (
                (
                    jac[i, j, k] * (
                        jac[i, j, k - 1] *
                        met[i, j, k - 1, 3, 3] *
                        turbulence_diffusion_coefficient(
                            state,
                            i,
                            j,
                            k - 1,
                            KM(),
                        )
                    ) +
                    jac[i, j, k - 1] * (
                        jac[i, j, k] *
                        met[i, j, k, 3, 3] *
                        turbulence_diffusion_coefficient(state, i, j, k, KM())
                    )
                ) / (jac[i, j, k] + jac[i, j, k - 1]) +
                (
                    jac[i, j + 1, k] * (
                        jac[i, j + 1, k - 1] *
                        met[i, j + 1, k - 1, 3, 3] *
                        turbulence_diffusion_coefficient(
                            state,
                            i,
                            j + 1,
                            k - 1,
                            KM(),
                        )
                    ) +
                    jac[i, j + 1, k - 1] * (
                        jac[i, j + 1, k] *
                        met[i, j + 1, k, 3, 3] *
                        turbulence_diffusion_coefficient(
                            state,
                            i,
                            j + 1,
                            k,
                            KM(),
                        )
                    )
                ) / (jac[i, j + 1, k] + jac[i, j + 1, k - 1])
            )

        jacc = 0.5 * (jac[i, j, k] + jac[i, j + 1, k])

        ith = i - nbx
        jth = j - nby
        kth = k - nbz

        ath[ith, jth, kth] = -dtdz2 * k33d / jacc
        bth[ith, jth, kth] = 1 + dtdz2 * k33u / jacc + dtdz2 * k33d / jacc
        cth[ith, jth, kth] = -dtdz2 * k33u / jacc
        fth[ith, jth, kth] =
            dtdz2 * k33u / jacc * v[i, j, k + 1] +
            (1 - dtdz2 * k33u / jacc - dtdz2 * k33d / jacc) * v[i, j, k] +
            dtdz2 * k33d / jacc * v[i, j, k - 1]
    end

    thomas_algorithm!(state)

    v[i0:i1, j0:j1, k0:k1] .= fth
    return
end

function turbulent_diffusion!(state::State, dt::AbstractFloat, variable::W)
    (; nbx, nby, nbz) = state.namelists.domain
    (; u, v, w) = state.variables.predictands
    (; uold, vold) = state.variables.backups
    (; i0, i1, j0, j1, k0, k1) = state.domain
    (; jac, met, dz) = state.grid
    (; ath, bth, cth, fth) = state.variables.auxiliaries

    dtdz2 = dt / (2.0 * dz^2.0)

    reset_thomas!(state)

    @ivy for k in k0:k1, j in j0:j1, i in i0:i1
        k33u =
            jac[i, j, k + 1] *
            met[i, j, k + 1, 3, 3] *
            turbulence_diffusion_coefficient(state, i, j, k + 1, KM())
        k33d =
            jac[i, j, k - 1] *
            met[i, j, k - 1, 3, 3] *
            turbulence_diffusion_coefficient(state, i, j, k - 1, KM())

        wu = compute_vertical_wind(i, j, k + 1, state)
        wc = compute_vertical_wind(i, j, k, state)
        wd = compute_vertical_wind(i, j, k - 1, state)

        dudtd =
            0.5 * (
                (u[i - 1, j, k] - uold[i - 1, j, k]) / dt +
                (u[i, j, k] - uold[i, j, k]) / dt
            )
        dudtu =
            0.5 * (
                (u[i - 1, j, k + 1] - uold[i - 1, j, k + 1]) / dt +
                (u[i, j, k + 1] - uold[i, j, k + 1]) / dt
            )
        dudtc13 =
            (
                jac[i, j, k] * met[i, j, k + 1, 1, 3] * dudtu +
                jac[i, j, k + 1] * met[i, j, k, 1, 3] * dudtd
            ) / (jac[i, j, k] + jac[i, j, k + 1])

        dvdtd =
            0.5 * (
                (v[i, j - 1, k] - vold[i, j - 1, k]) / dt +
                (v[i, j, k] - vold[i, j, k]) / dt
            )
        dvdtu =
            0.5 * (
                (v[i, j - 1, k + 1] - vold[i, j - 1, k + 1]) / dt +
                (v[i, j, k + 1] - vold[i, j, k + 1]) / dt
            )
        dvdtc23 =
            (
                jac[i, j, k] * met[i, j, k + 1, 2, 3] * dvdtu +
                jac[i, j, k + 1] * met[i, j, k, 2, 3] * dvdtd
            ) / (jac[i, j, k] + jac[i, j, k + 1])

        jacc =
            2.0 * jac[i, j, k] * jac[i, j, k + 1] /
            (jac[i, j, k] + jac[i, j, k + 1])

        ith = i - nbx
        jth = j - nby
        kth = k - nbz

        ath[ith, jth, kth] = -dtdz2 / jacc * k33d
        bth[ith, jth, kth] = 1 + dtdz2 / jacc * k33u + dtdz2 / jacc * k33d
        cth[ith, jth, kth] = -dtdz2 / jacc * k33u
        fth[ith, jth, kth] =
            dtdz2 / jacc * k33u * wu +
            (1 - dtdz2 / jacc * k33u - dtdz2 / jacc * k33d) * wc +
            dtdz2 / jacc * k33d * wd +
            dudtc13 +
            dvdtc23
    end

    thomas_algorithm!(state)

    w[i0:i1, j0:j1, k0:k1] .= fth
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
    (; p) = state.variables.predictands
    (; i0, i1, j0, j1, k0, k1) = state.domain
    (; nbx, nby, nbz) = state.namelists.domain
    (; jac, met, dz) = state.grid
    (; ath, bth, cth, fth) = state.variables.auxiliaries

    dtdz2 = dt / (2.0 * dz^2.0)

    reset_thomas!(state)

    @ivy for k in k0:k1, j in j0:j1, i in i0:i1
        kh33d =
            (
                jac[i, j, k - 1] * (
                    jac[i, j, k] *
                    met[i, j, k, 3, 3] *
                    turbulence_diffusion_coefficient(state, i, j, k, KH())
                ) +
                jac[i, j, k] * (
                    jac[i, j, k - 1] *
                    met[i, j, k - 1, 3, 3] *
                    turbulence_diffusion_coefficient(state, i, j, k - 1, KH())
                )
            ) / (jac[i, j, k - 1] + jac[i, j, k])
        kh33u =
            (
                jac[i, j, k + 1] * (
                    jac[i, j, k] *
                    met[i, j, k] *
                    turbulence_diffusion_coefficient(state, i, j, k, KH())
                ) +
                jac[i, j, k] * (
                    jac[i, j, k + 1] *
                    met[i, j, k + 1, 3, 3] *
                    turbulence_diffusion_coefficient(state, i, j, k + 1, KH())
                )
            ) / (jac[i, j, k + 1] + jac[i, j, k])

        ith = i - nbx
        jth = j - nby
        kth = k - nbz

        ath[ith, jth, kth] = -dtdz2 / jac[i, j, k] * kh33d
        bth[ith, jth, kth] =
            1 + dtdz2 / jac[i, j, k] * kh33u + dtdz2 / jac[i, j, k] * kh33d
        cth[ith, jth, kth] = -dtdz2 / jac[i, j, k] * kh33u

        fth[ith, jth, kth] =
            (1 - dtdz2 / jac[i, j, k] * kh33u - dtdz2 / jac[i, j, k] * kh33d) *
            p[i, j, k] +
            dtdz2 / jac[i, j, k] * kh33u * p[i, j, k + 1] +
            dtdz2 / jac[i, j, k] * kh33d * p[i, j, k - 1]
    end

    thomas_algorithm!(state)

    p[i0:i1, j0:j1, k0:k1] .= fth
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
    (; i0, i1, j0, j1, k0, k1) = state.domain
    (; nbx, nby, nbz) = state.namelists.domain
    (; jac, met, dz) = state.grid
    (; ath, bth, cth, fth) = state.variables.auxiliaries

    dtdz2 = dt / (2.0 * dz^2.0)

    reset_thomas!(state)

    for field in 1:fieldcount(TracerPredictands)
        chi = getfield(tracerpredictands, field)
        @ivy for k in k0:k1, j in j0:j1, i in i0:i1
            kh33d =
                (
                    jac[i, j, k - 1] * (
                        jac[i, j, k] *
                        met[i, j, k, 3, 3] *
                        turbulence_diffusion_coefficient(state, i, j, k, KH())
                    ) +
                    jac[i, j, k] * (
                        jac[i, j, k - 1] *
                        met[i, j, k - 1, 3, 3] *
                        turbulence_diffusion_coefficient(
                            state,
                            i,
                            j,
                            k - 1,
                            KH(),
                        )
                    )
                ) / (jac[i, j, k - 1] + jac[i, j, k])
            kh33u =
                (
                    jac[i, j, k + 1] * (
                        jac[i, j, k] *
                        met[i, j, k] *
                        turbulence_diffusion_coefficient(state, i, j, k, KH())
                    ) +
                    jac[i, j, k] * (
                        jac[i, j, k + 1] *
                        met[i, j, k + 1, 3, 3] *
                        turbulence_diffusion_coefficient(
                            state,
                            i,
                            j,
                            k + 1,
                            KH(),
                        )
                    )
                ) / (jac[i, j, k + 1] + jac[i, j, k])

            ith = i - nbx
            jth = j - nby
            kth = k - nbz

            ath[ith, jth, kth] = -dtdz2 / jac[i, j, k] * kh33d
            bth[ith, jth, kth] =
                1 + dtdz2 / jac[i, j, k] * kh33u + dtdz2 / jac[i, j, k] * kh33d
            cth[ith, jth, kth] = -dtdz2 / jac[i, j, k] * kh33u

            fth[ith, jth, kth] =
                (
                    1 - dtdz2 / jac[i, j, k] * kh33u -
                    dtdz2 / jac[i, j, k] * kh33d
                ) * chi[i, j, k] +
                dtdz2 / jac[i, j, k] * kh33u * chi[i, j, k + 1] +
                dtdz2 / jac[i, j, k] * kh33d * chi[i, j, k - 1]
        end

        thomas_algorithm!(state)

        chi[i0:i1, j0:j1, k0:k1] .= fth
    end
    return
end
