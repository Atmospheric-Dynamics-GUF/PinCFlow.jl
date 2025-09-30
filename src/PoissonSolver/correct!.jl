"""
```julia
correct!(state::State, dt::AbstractFloat, rayleigh_factor::AbstractFloat)
```

Correct the Exner-pressure, wind and density fluctuations such that the divergence constraint is satisfied, using the Exner-pressure differences obtained from the solution to the Poisson problem.

This method calls specific methods for the correction of each individual variable.

```julia
correct!(
    state::State,
    dt::AbstractFloat,
    variable::U,
    rayleigh_factor::AbstractFloat,
)
```

Correct the zonal wind to account for the pressure differences obtained from the solution to the Poisson problem.

The correction is given by

```math
u_{i + 1 / 2} \\rightarrow u_{i + 1 / 2} - \\mathcal{C}_{i + 1 / 2}^{\\rho u}
```

in Boussinesq/pseudo-incompressible mode and

```math
U_{i + 1 / 2} \\rightarrow U_{i + 1 / 2} - \\left(J P\\right)_{i + 1 / 2} \\mathcal{C}_{i + 1 / 2}^{\\rho u}
```

in compressible mode, with

```math
\\mathcal{C}_{i + 1 / 2}^{\\rho u} = \\left(1 + \\beta_{\\mathrm{R}, i + 1 / 2} \\Delta t\\right)^{- 1} \\Delta t c_p \\frac{P_{i + 1 / 2}}{\\rho_{i + 1 / 2}} \\mathcal{D}_{i + 1 / 2}^{\\rho u}.
```

Therein, ``\\Delta t`` is the fractional time step given as input to this method and ``c_p \\left(P_{i + 1 / 2} / \\rho_{i + 1 / 2}\\right) \\mathcal{D}_{i + 1 / 2}^{\\rho u}`` is computed with `compute_pressure_gradient`.

```julia
correct!(
    state::State,
    dt::AbstractFloat,
    variable::V,
    rayleigh_factor::AbstractFloat,
)
```

Correct the meridional wind to account for the pressure differences obtained from the solution to the Poisson problem.

The correction is given by

```math
v_{j + 1 / 2} \\rightarrow v_{j + 1 / 2} - \\mathcal{C}_{j + 1 / 2}^{\\rho v}
```

in Boussinesq/pseudo-incompressible mode and

```math
V_{j + 1 / 2} \\rightarrow V_{j + 1 / 2} - \\left(J P\\right)_{j + 1 / 2} \\mathcal{C}_{j + 1 / 2}^{\\rho v}
```

in compressible mode, with

```math
\\mathcal{C}_{j + 1 / 2}^{\\rho v} = \\left(1 + \\beta_{\\mathrm{R}, j + 1 / 2} \\Delta t\\right)^{- 1} \\Delta t c_p \\frac{P_{j + 1 / 2}}{\\rho_{j + 1 / 2}} \\mathcal{D}_{j + 1 / 2}^{\\rho v},
```

where ``c_p \\left(P_{j + 1 / 2} / \\rho_{j + 1 / 2}\\right) \\mathcal{D}_{j + 1 / 2}^{\\rho v}`` is computed with `compute_pressure_gradient`.

```julia
correct!(
    state::State,
    dt::AbstractFloat,
    variable::W,
    rayleigh_factor::AbstractFloat,
)
```

Correct the transformed vertical wind to account for the pressure differences obtained from the solution to the Poisson problem.

The correction is given by

```math
\\begin{align*}
    \\widehat{w}_{k + 1 / 2} & \\rightarrow \\widehat{w}_{k + 1 / 2} - \\left[1 + \\beta_{\\mathrm{R}, k + 1 / 2} \\Delta t + \\frac{\\overline{\\rho}_{k + 1 / 2}}{\\rho_{k + 1 / 2}} \\left(N \\Delta t\\right)^2\\right]^{- 1}\\\\
    & \\quad \\times \\left\\{\\Delta t c_p \\frac{P_{k + 1 / 2}}{\\rho_{k + 1 / 2}} \\mathcal{D}_{k + 1 / 2}^{\\rho \\widehat{w}} + \\frac{\\overline{\\rho}_{k + 1 / 2}}{\\rho_{k + 1 / 2}} \\left(N \\Delta t\\right)^2 \\left[\\left(G^{1 3} \\mathcal{C}^{\\rho u}\\right)_{k + 1 / 2} + \\left(G^{23} \\mathcal{C}^{\\rho v}\\right)_{k + 1 / 2}\\right]\\right\\},
\\end{align*}
```

in Boussinesq/pseudo-incompressible mode and

```math
\\begin{align*}
    \\widehat{W}_{k + 1 / 2} & \\rightarrow \\widehat{W}_{k + 1 / 2} - \\left[1 + \\beta_{\\mathrm{R}, k + 1 / 2} \\Delta t + \\frac{\\left(P / \\overline{\\theta}\\right)_{k + 1 / 2}}{\\rho_{k + 1 / 2}} \\left(N \\Delta t\\right)^2\\right]^{- 1}\\\\
    & \\quad \\times \\left\\{\\Delta t c_p \\left(J P\\right)_{k + 1 / 2} \\frac{P_{k + 1 / 2}}{\\rho_{k + 1 / 2}} \\mathcal{D}_{k + 1 / 2}^{\\rho \\widehat{w}} \\vphantom{\\frac{\\left(P / \\overline{\\theta}\\right)_{k + 1 / 2}}{\\rho_{k + 1 / 2}}}\\right.\\\\
    & \\qquad \\quad + \\left.\\left(J P\\right)_{k + 1 / 2} \\frac{\\left(P / \\overline{\\theta}\\right)_{k + 1 / 2}}{\\rho_{k + 1 / 2}} \\left(N \\Delta t\\right)^2 \\left[\\left(G^{1 3} \\mathcal{C}^{\\rho u}\\right)_{k + 1 / 2} + \\left(G^{23} \\mathcal{C}^{\\rho v}\\right)_{k + 1 / 2}\\right]\\right\\},
\\end{align*}
```

in compressible mode, where ``c_p \\left(P_{k + 1 / 2} / \\rho_{k + 1 / 2}\\right) \\mathcal{D}_{k + 1 / 2}^{\\rho \\widehat{w}}`` is computed with `compute_pressure_gradient`.

```julia
correct!(
    state::State,
    dt::AbstractFloat,
    variable::RhoP,
    rayleigh_factor::AbstractFloat,
)
```

Correct the density fluctuations to account for the pressure differences obtained from the solution to the Poisson problem.

The correction is given by

```math
\\begin{align*}
    \\rho' & \\rightarrow \\rho' + \\frac{\\rho}{g} \\left[1 + \\beta_\\mathrm{R} \\Delta t + \\frac{\\overline{\\rho}}{\\rho} \\left(N \\Delta t\\right)^2\\right]^{- 1}\\\\
    & \\quad \\times \\left[- \\frac{\\overline{\\rho}}{\\rho} \\left(N \\Delta t\\right)^2 J \\left(c_p \\frac{P_{k + 1 / 2}}{\\rho_{k + 1 / 2}} \\mathcal{D}_{k + 1 / 2}^{\\rho \\widehat{w}}\\right)\\right.\\\\
    & \\qquad \\quad + \\left.\\frac{\\overline{\\rho}}{\\rho} N^2 \\Delta t J \\left(1 + \\beta_\\mathrm{R} \\Delta t\\right) \\left(G^{1 3} \\mathcal{C}^{\\rho u} + G^{2 3} \\mathcal{C}^{\\rho v}\\right)\\right],
\\end{align*}
```

in Boussinesq/pseudo-incompressible mode and

```math
\\begin{align*}
    \\rho' & \\rightarrow \\rho' + \\frac{\\rho}{g} \\left[1 + \\beta_\\mathrm{R} \\Delta t + \\frac{P / \\overline{\\theta}}{\\rho} \\left(N \\Delta t\\right)^2\\right]^{- 1}\\\\
    & \\quad \\times \\left[- \\frac{P / \\overline{\\theta}}{\\rho} \\left(N \\Delta t\\right)^2 J \\left(c_p \\frac{P_{k + 1 / 2}}{\\rho_{k + 1 / 2}} \\mathcal{D}_{k + 1 / 2}^{\\rho \\widehat{w}}\\right)\\right.\\\\
    & \\qquad \\quad + \\left.\\frac{P / \\overline{\\theta}}{\\rho} N^2 \\Delta t J \\left(1 + \\beta_\\mathrm{R} \\Delta t\\right) \\left(G^{1 3} \\mathcal{C}^{\\rho u} + G^{2 3} \\mathcal{C}^{\\rho v}\\right)\\right],
\\end{align*}
```

in compressible mode, where ``c_p \\left(P_{k + 1 / 2} / \\rho_{k + 1 / 2}\\right) \\mathcal{D}_{k + 1 / 2}^{\\rho \\widehat{w}}`` and ``c_p \\left(P_{k - 1 / 2} / \\rho_{k - 1 / 2}\\right) \\mathcal{D}_{k - 1 / 2}^{\\rho \\widehat{w}}`` are computed with `compute_pressure_gradient`, and used to interpolate to ``\\left(i, j, k\\right)``.

```julia
correct!(state::State, variable::PiP)
```

Update the Exner-pressure fluctuations with the differences obtained from the solution to the Poisson problem.

# Arguments

  - `state`: Model state.

  - `dt`: Time step.

  - `variable`: Variable to correct.

  - `rayleigh_factor`: Factor by which the Rayleigh-damping coefficient is multiplied.

# See also

  - [`PinCFlow.Update.compute_pressure_gradient`](@ref)

  - [`PinCFlow.Update.compute_compressible_wind_factor`](@ref)

  - [`PinCFlow.Update.compute_buoyancy_factor`](@ref)
"""
function correct! end

function correct!(
    state::State,
    dt::AbstractFloat,
    rayleigh_factor::AbstractFloat,
)
    correct!(state, dt, U(), rayleigh_factor)
    correct!(state, dt, V(), rayleigh_factor)
    correct!(state, dt, W(), rayleigh_factor)
    correct!(state, dt, RhoP(), rayleigh_factor)
    correct!(state, PiP())
    return
end

function correct!(
    state::State,
    dt::AbstractFloat,
    variable::U,
    rayleigh_factor::AbstractFloat,
)
    (; spongelayer, sponge_uv) = state.namelists.sponge
    (; ndzz, nzz, ko, i0, i1, j0, j1, k0, k1) = state.domain
    (; betar) = state.sponge
    (; corx) = state.poisson.correction
    (; dpip) = state.variables.increments
    (; u) = state.variables.predictands

    kmin = k0
    kmax = ko + nzz == ndzz ? k1 : k1 + 1

    @ivy for k in kmin:kmax, j in j0:j1, i in (i0 - 1):i1
        factor = 1.0

        if spongelayer && sponge_uv
            factor +=
                dt *
                0.5 *
                (betar[i, j, k] + betar[i + 1, j, k]) *
                rayleigh_factor
        end

        gradient = compute_pressure_gradient(state, dpip, i, j, k, U())

        corx[i, j, k] = dt / factor * gradient

        jpedger = compute_compressible_wind_factor(state, i, j, k, U())

        u[i, j, k] -= jpedger * corx[i, j, k]
    end

    return
end

function correct!(
    state::State,
    dt::AbstractFloat,
    variable::V,
    rayleigh_factor::AbstractFloat,
)
    (; spongelayer, sponge_uv) = state.namelists.sponge
    (; ndzz, nzz, ko, i0, i1, j0, j1, k0, k1) = state.domain
    (; betar) = state.sponge
    (; cory) = state.poisson.correction
    (; dpip) = state.variables.increments
    (; v) = state.variables.predictands

    kmin = k0
    kmax = ko + nzz == ndzz ? k1 : k1 + 1

    @ivy for k in kmin:kmax, j in (j0 - 1):j1, i in i0:i1
        factor = 1.0

        if spongelayer && sponge_uv
            factor +=
                dt *
                0.5 *
                (betar[i, j, k] + betar[i, j + 1, k]) *
                rayleigh_factor
        end

        gradient = compute_pressure_gradient(state, dpip, i, j, k, V())

        cory[i, j, k] = dt / factor * gradient

        jpedgef = compute_compressible_wind_factor(state, i, j, k, V())

        v[i, j, k] -= jpedgef * cory[i, j, k]
    end

    return
end

function correct!(
    state::State,
    dt::AbstractFloat,
    variable::W,
    rayleigh_factor::AbstractFloat,
)
    (; spongelayer) = state.namelists.sponge
    (; ndzz, nzz, ko, i0, i1, j0, j1, k0, k1) = state.domain
    (; jac, met) = state.grid
    (; bvsstrattfc) = state.atmosphere
    (; betar) = state.sponge
    (; corx, cory) = state.poisson.correction
    (; dpip) = state.variables.increments
    (; w) = state.variables.predictands

    kmin = ko == 0 ? k0 : k0 - 1
    kmax = ko + nzz == ndzz ? k1 - 1 : k1

    @ivy for k in kmin:kmax, j in j0:j1, i in i0:i1
        factor = 1.0

        if spongelayer
            factor +=
                dt * (
                    jac[i, j, k + 1] * betar[i, j, k] +
                    jac[i, j, k] * betar[i, j, k + 1]
                ) / (jac[i, j, k] + jac[i, j, k + 1]) * rayleigh_factor
        end

        bvsstratedgeu =
            (
                jac[i, j, k + 1] * bvsstrattfc[i, j, k] +
                jac[i, j, k] * bvsstrattfc[i, j, k + 1]
            ) / (jac[i, j, k] + jac[i, j, k + 1])

        gradient = compute_pressure_gradient(state, dpip, i, j, k, W())

        jpedgeu = compute_compressible_wind_factor(state, i, j, k, W())
        fw = compute_buoyancy_factor(state, i, j, k, W())

        w[i, j, k] +=
            -dt / (factor + fw * bvsstratedgeu * dt^2.0) * jpedgeu * gradient -
            1.0 / (factor + fw * bvsstratedgeu * dt^2.0) *
            fw *
            bvsstratedgeu *
            dt^2.0 *
            jpedgeu *
            0.5 *
            (
                jac[i, j, k + 1] * (
                    met[i, j, k, 1, 3] * (corx[i, j, k] + corx[i - 1, j, k]) +
                    met[i, j, k, 2, 3] * (cory[i, j, k] + cory[i, j - 1, k])
                ) +
                jac[i, j, k] * (
                    met[i, j, k + 1, 1, 3] *
                    (corx[i, j, k + 1] + corx[i - 1, j, k + 1]) +
                    met[i, j, k + 1, 2, 3] *
                    (cory[i, j, k + 1] + cory[i, j - 1, k + 1])
                )
            ) / (jac[i, j, k] + jac[i, j, k + 1])
    end

    return
end

function correct!(
    state::State,
    dt::AbstractFloat,
    variable::RhoP,
    rayleigh_factor::AbstractFloat,
)
    (; nbz) = state.namelists.domain
    (; spongelayer) = state.namelists.sponge
    (; g_ndim) = state.constants
    (; ndzz, ko, i0, i1, j0, j1, k0, k1) = state.domain
    (; jac, met) = state.grid
    (; rhostrattfc, bvsstrattfc) = state.atmosphere
    (; betar) = state.sponge
    (; corx, cory) = state.poisson.correction
    (; dpip) = state.variables.increments
    (; rho, rhop) = state.variables.predictands

    @ivy for k in k0:k1, j in j0:j1, i in i0:i1
        factor = 1.0

        if spongelayer
            factor += dt * betar[i, j, k] * rayleigh_factor
        end

        lower_gradient =
            compute_pressure_gradient(state, dpip, i, j, k - 1, W())
        upper_gradient = compute_pressure_gradient(state, dpip, i, j, k, W())

        if ko + k == k0
            lower_gradient = 0.0
        elseif ko + k == ndzz - nbz
            upper_gradient = 0.0
        end

        gradient = 0.5 * (lower_gradient + upper_gradient)

        fb = compute_buoyancy_factor(state, i, j, k, RhoP())
        db =
            -1.0 / (factor + fb * bvsstrattfc[i, j, k] * dt^2.0) * (
                -fb * bvsstrattfc[i, j, k] * dt^2.0 * jac[i, j, k] * gradient +
                fb *
                bvsstrattfc[i, j, k] *
                dt *
                jac[i, j, k] *
                factor *
                0.5 *
                (
                    met[i, j, k, 1, 3] * (corx[i, j, k] + corx[i - 1, j, k]) +
                    met[i, j, k, 2, 3] * (cory[i, j, k] + cory[i, j - 1, k])
                )
            )

        rhop[i, j, k] -= (rho[i, j, k] + rhostrattfc[i, j, k]) / g_ndim * db
    end

    return
end

function correct!(state::State, variable::PiP)
    (; i0, i1, j0, j1, k0, k1) = state.domain
    (; pip) = state.variables.predictands
    (; dpip) = state.variables.increments

    ii = (i0 - 1):(i1 + 1)
    jj = (j0 - 1):(j1 + 1)
    kk = (k0 - 1):(k1 + 1)

    pip[ii, jj, kk] .+= dpip[ii, jj, kk]

    return
end
