"""
```julia
compute_pressure_gradient(
    state::State,
    pip::AbstractArray{<:AbstractFloat, 3},
    i::Integer,
    j::Integer,
    k::Integer,
    variable::U,
)::AbstractFloat
```

Compute and return the pressure(-difference)-gradient term in the zonal-wind equation at ``\\left(i + 1 / 2, j, k\\right)``, using the pressure(-difference) field `pip`.

The pressure-gradient component is given by

```math
\\mathcal{P}^{\\rho u}_{i + 1 / 2} = \\frac{\\pi'_{i + 1} - \\pi'}{\\Delta \\widehat{x}} + G^{13}_{i + 1 / 2} \\frac{\\pi'_{i + 1 / 2, k + 1} - \\pi'_{i + 1 / 2, k - 1}}{2 \\Delta \\widehat{z}}.
```

Since the Exner-pressure is not known in the vertical ghost cells, a different discretization is needed at the vertical boundaries. At ``k = k_0`` (in the first process in ``\\widehat{z}``), the alternative second-order-accurate approximation

```math
\\mathcal{P}^{\\rho u}_{i + 1 / 2} = \\frac{\\pi'_{i + 1} - \\pi'}{\\Delta \\widehat{x}} + G^{13}_{i + 1 / 2} \\frac{- \\pi'_{i + 1 / 2, k + 2} + 4 \\pi'_{i + 1 / 2, k + 1} - 3 \\pi'_{i + 1 / 2}}{2 \\Delta \\widehat{z}}
```

is used and, in a similar manner, one has

```math
\\mathcal{P}^{\\rho u}_{i + 1 / 2} = \\frac{\\pi'_{i + 1} - \\pi'}{\\Delta \\widehat{x}} + G^{13}_{i + 1 / 2} \\frac{\\pi'_{i + 1 / 2, k - 2} - 4 \\pi'_{i + 1 / 2, k - 1} + 3 \\pi'_{i + 1 / 2}}{2 \\Delta \\widehat{z}}
```

at ``k = k_1`` (in the last process in ``\\widehat{z}``). The corresponding pressure-difference-gradient component ``\\mathcal{D}^{\\rho u}_{i + 1 / 2}`` is obtained by replacing ``\\pi'`` with ``\\Delta \\pi'``. The returned quantity also includes the factor ``c_p \\left(P_{i + 1 / 2} / \\rho_{i + 1 / 2}\\right)``.

```julia
compute_pressure_gradient(
    state::State,
    pip::AbstractArray{<:AbstractFloat, 3},
    i::Integer,
    j::Integer,
    k::Integer,
    variable::V,
)::AbstractFloat
```

Compute and return the pressure-gradient term in the meridional-wind equation at ``\\left(i, j + 1 / 2, k\\right)``, using the pressure(-difference) field `pip`.

The pressure-gradient component is given by

```math
\\mathcal{P}^{\\rho v}_{j + 1 / 2} = \\frac{\\pi'_{j + 1} - \\pi'}{\\Delta \\widehat{y}} + G^{23}_{j + 1 / 2} \\frac{\\pi'_{j + 1 / 2, k + 1} - \\pi'_{j + 1 / 2, k - 1}}{2 \\Delta \\widehat{z}}.
```

Analogous to the component in the zonal-wind equation, one has

```math
\\mathcal{P}^{\\rho v}_{j + 1 / 2} = \\frac{\\pi'_{j + 1} - \\pi'}{\\Delta \\widehat{y}} + G^{23}_{j + 1 / 2} \\frac{- \\pi'_{j + 1 / 2, k + 2} + 4 \\pi'_{j + 1 / 2, k + 1} - 3 \\pi'_{j + 1 / 2}}{2 \\Delta \\widehat{z}}
```

at ``k = k_0`` (in the first process in ``\\widehat{z}``) and

```math
\\mathcal{P}^{\\rho v}_{j + 1 / 2} = \\frac{\\pi'_{j + 1} - \\pi'}{\\Delta \\widehat{y}} + G^{23}_{j + 1 / 2} \\frac{\\pi'_{j + 1 / 2, k - 2} - 4 \\pi'_{j + 1 / 2, k - 1} + 3 \\pi'_{j + 1 / 2}}{2 \\Delta \\widehat{z}}
```

at ``k = k_1`` (in the last process in ``\\widehat{z}``). The corresponding pressure-difference-gradient component ``\\mathcal{D}^{\\rho v}_{j + 1 / 2}`` is obtained by replacing ``\\pi'`` with ``\\Delta \\pi'``. The returned quantity also includes the factor ``c_p \\left(P_{j + 1 / 2} / \\rho_{j + 1 / 2}\\right)``.

```julia
compute_pressure_gradient(
    state::State,
    pip::AbstractArray{<:AbstractFloat, 3},
    i::Integer,
    j::Integer,
    k::Integer,
    variable::W,
)::AbstractFloat
```

Compute and return the pressure(-difference)-gradient term in the transformed-vertical-wind equation at ``\\left(i, j, k + 1 / 2\\right)``, using the pressure(-difference) field `pip`.

The pressure-gradient component is given by

```math
\\begin{align*}
    \\mathcal{P}^{\\rho \\widehat{w}}_{k + 1 / 2} & = G^{13}_{k + 1 / 2} \\frac{\\pi'_{i + 1, k + 1 / 2} - \\pi'_{i - 1, k + 1 / 2}}{2 \\Delta \\widehat{x}} + G^{23}_{k + 1 / 2} \\frac{\\pi'_{j + 1, k + 1 / 2} - \\pi'_{j - 1, k + 1 / 2}}{2 \\Delta \\widehat{y}}\\\\
    & \\quad + G^{33}_{k + 1 / 2} \\frac{\\pi'_{k + 1} - \\pi'}{\\Delta \\widehat{z}}.
 \\end{align*}
```

At ``k = k_0 - 1`` (in the first process in ``\\widehat{z}``) and ``k = k_1`` (in the last process in ``\\widehat{z}``), it is set to zero. The corresponding pressure-difference-gradient component ``\\mathcal{D}^{\\rho \\widehat{w}}_{k + 1 / 2}`` is obtained by replacing ``\\pi'`` with ``\\Delta \\pi'``. The returned quantity also includes the factor ``c_p \\left(P_{k + 1 / 2} / \\rho_{k + 1 / 2}\\right)``.

# Arguments

  - `state`: Model state.

  - `pip`: Pressure field.

  - `i`: Zonal grid-cell index.

  - `j`: Meridional grid-cell index.

  - `k`: Vertical grid-cell index.

  - `variable`: Equation in which the respective pressure-gradient component is needed.
"""
function compute_pressure_gradient end

function compute_pressure_gradient(
    state::State,
    pip::AbstractArray{<:AbstractFloat, 3},
    i::Integer,
    j::Integer,
    k::Integer,
    variable::U,
)::AbstractFloat
    (; nbz) = state.namelists.domain
    (; kappainv, mainv2) = state.constants
    (; zz_size, ko, k0) = state.domain
    (; dx, dz, met) = state.grid
    (; rhostrattfc, pstrattfc) = state.atmosphere
    (; rho) = state.variables.predictands

    # Interpolate the density, mass-weighted potential temperature and metric
    # tensor element.
    @ivy rhoedger = 0.5 * (rho[i, j, k] + rho[i + 1, j, k])
    @ivy rhostratedger = 0.5 * (rhostrattfc[i, j, k] + rhostrattfc[i + 1, j, k])
    @ivy rhoedger += rhostratedger
    @ivy pedger = 0.5 * (pstrattfc[i, j, k] + pstrattfc[i + 1, j, k])
    @ivy met13edger = 0.5 * (met[i, j, k, 1, 3] + met[i + 1, j, k, 1, 3])

    # Compute the pressure gradient component.
    @ivy if ko + k == k0
        pipuuedger = 0.5 * (pip[i, j, k + 2] + pip[i + 1, j, k + 2])
        pipuedger = 0.5 * (pip[i, j, k + 1] + pip[i + 1, j, k + 1])
        pipedger = 0.5 * (pip[i, j, k] + pip[i + 1, j, k])
        gradient =
            kappainv * mainv2 * pedger / rhoedger * (
                (pip[i + 1, j, k] - pip[i, j, k]) / dx +
                met13edger *
                (-pipuuedger + 4.0 * pipuedger - 3.0 * pipedger) *
                0.5 / dz
            )
    elseif ko + k == zz_size - nbz
        pipddedger = 0.5 * (pip[i, j, k - 2] + pip[i + 1, j, k - 2])
        pipdedger = 0.5 * (pip[i, j, k - 1] + pip[i + 1, j, k - 1])
        pipedger = 0.5 * (pip[i, j, k] + pip[i + 1, j, k])
        gradient =
            kappainv * mainv2 * pedger / rhoedger * (
                (pip[i + 1, j, k] - pip[i, j, k]) / dx +
                met13edger *
                (pipddedger - 4.0 * pipdedger + 3.0 * pipedger) *
                0.5 / dz
            )
    else
        pipuedger = 0.5 * (pip[i, j, k + 1] + pip[i + 1, j, k + 1])
        pipdedger = 0.5 * (pip[i, j, k - 1] + pip[i + 1, j, k - 1])
        gradient =
            kappainv * mainv2 * pedger / rhoedger * (
                (pip[i + 1, j, k] - pip[i, j, k]) / dx +
                met13edger * (pipuedger - pipdedger) * 0.5 / dz
            )
    end

    return gradient
end

function compute_pressure_gradient(
    state::State,
    pip::AbstractArray{<:AbstractFloat, 3},
    i::Integer,
    j::Integer,
    k::Integer,
    variable::V,
)::AbstractFloat
    (; nbz) = state.namelists.domain
    (; kappainv, mainv2) = state.constants
    (; zz_size, ko, k0) = state.domain
    (; dy, dz, met) = state.grid
    (; rhostrattfc, pstrattfc) = state.atmosphere
    (; rho) = state.variables.predictands

    # Interpolate the density, mass-weighted potential temperature and metric
    # tensor element.
    @ivy rhoedgef = 0.5 * (rho[i, j, k] + rho[i, j + 1, k])
    @ivy rhostratedgef = 0.5 * (rhostrattfc[i, j, k] + rhostrattfc[i, j + 1, k])
    @ivy rhoedgef += rhostratedgef
    @ivy pedgef = 0.5 * (pstrattfc[i, j, k] + pstrattfc[i, j + 1, k])
    @ivy met23edgef = 0.5 * (met[i, j, k, 2, 3] + met[i, j + 1, k, 2, 3])

    # Compute the pressure gradient component.
    @ivy if ko + k == k0
        pipuuedgef = 0.5 * (pip[i, j, k + 2] + pip[i, j + 1, k + 2])
        pipuedgef = 0.5 * (pip[i, j, k + 1] + pip[i, j + 1, k + 1])
        pipedgef = 0.5 * (pip[i, j, k] + pip[i, j + 1, k])
        gradient =
            kappainv * mainv2 * pedgef / rhoedgef * (
                (pip[i, j + 1, k] - pip[i, j, k]) / dy +
                met23edgef *
                (-pipuuedgef + 4.0 * pipuedgef - 3.0 * pipedgef) *
                0.5 / dz
            )
    elseif ko + k == zz_size - nbz
        pipddedgef = 0.5 * (pip[i, j, k - 2] + pip[i, j + 1, k - 2])
        pipdedgef = 0.5 * (pip[i, j, k - 1] + pip[i, j + 1, k - 1])
        pipedgef = 0.5 * (pip[i, j, k] + pip[i, j + 1, k])
        gradient =
            kappainv * mainv2 * pedgef / rhoedgef * (
                (pip[i, j + 1, k] - pip[i, j, k]) / dy +
                met23edgef *
                (pipddedgef - 4.0 * pipdedgef + 3.0 * pipedgef) *
                0.5 / dz
            )
    else
        pipuedgef = 0.5 * (pip[i, j, k + 1] + pip[i, j + 1, k + 1])
        pipdedgef = 0.5 * (pip[i, j, k - 1] + pip[i, j + 1, k - 1])
        gradient =
            kappainv * mainv2 * pedgef / rhoedgef * (
                (pip[i, j + 1, k] - pip[i, j, k]) / dy +
                met23edgef * (pipuedgef - pipdedgef) * 0.5 / dz
            )
    end

    return gradient
end

function compute_pressure_gradient(
    state::State,
    pip::AbstractArray{<:AbstractFloat, 3},
    i::Integer,
    j::Integer,
    k::Integer,
    variable::W,
)::AbstractFloat
    (; kappainv, mainv2) = state.constants
    (; dx, dy, dz, jac, met) = state.grid
    (; rhostrattfc, pstrattfc) = state.atmosphere
    (; rho) = state.variables.predictands

    # Interpolate the density, mass-weighted potential temperature and metric
    # tensor element.
    @ivy rhoedgeu =
        (jac[i, j, k + 1] * rho[i, j, k] + jac[i, j, k] * rho[i, j, k + 1]) /
        (jac[i, j, k] + jac[i, j, k + 1])
    @ivy rhoedgeu +=
        (
            jac[i, j, k + 1] * rhostrattfc[i, j, k] +
            jac[i, j, k] * rhostrattfc[i, j, k + 1]
        ) / (jac[i, j, k] + jac[i, j, k + 1])
    @ivy pedgeu =
        (
            jac[i, j, k + 1] * pstrattfc[i, j, k] +
            jac[i, j, k] * pstrattfc[i, j, k + 1]
        ) / (jac[i, j, k] + jac[i, j, k + 1])
    @ivy met13edgeu =
        (
            jac[i, j, k + 1] * met[i, j, k, 1, 3] +
            jac[i, j, k] * met[i, j, k + 1, 1, 3]
        ) / (jac[i, j, k] + jac[i, j, k + 1])
    @ivy met23edgeu =
        (
            jac[i, j, k + 1] * met[i, j, k, 2, 3] +
            jac[i, j, k] * met[i, j, k + 1, 2, 3]
        ) / (jac[i, j, k] + jac[i, j, k + 1])
    @ivy met33edgeu =
        (
            jac[i, j, k + 1] * met[i, j, k, 3, 3] +
            jac[i, j, k] * met[i, j, k + 1, 3, 3]
        ) / (jac[i, j, k] + jac[i, j, k + 1])

    # Compute the pressure gradient component.
    @ivy pipredgeu =
        (
            jac[i + 1, j, k + 1] * pip[i + 1, j, k] +
            jac[i + 1, j, k] * pip[i + 1, j, k + 1]
        ) / (jac[i + 1, j, k] + jac[i + 1, j, k + 1])
    @ivy pipledgeu =
        (
            jac[i - 1, j, k + 1] * pip[i - 1, j, k] +
            jac[i - 1, j, k] * pip[i - 1, j, k + 1]
        ) / (jac[i - 1, j, k] + jac[i - 1, j, k + 1])
    @ivy pipfedgeu =
        (
            jac[i, j + 1, k + 1] * pip[i, j + 1, k] +
            jac[i, j + 1, k] * pip[i, j + 1, k + 1]
        ) / (jac[i, j + 1, k] + jac[i, j + 1, k + 1])
    @ivy pipbedgeu =
        (
            jac[i, j - 1, k + 1] * pip[i, j - 1, k] +
            jac[i, j - 1, k] * pip[i, j - 1, k + 1]
        ) / (jac[i, j - 1, k] + jac[i, j - 1, k + 1])
    @ivy gradient =
        kappainv * mainv2 * pedgeu / rhoedgeu * (
            met13edgeu * (pipredgeu - pipledgeu) * 0.5 / dx +
            met23edgeu * (pipfedgeu - pipbedgeu) * 0.5 / dy +
            met33edgeu * (pip[i, j, k + 1] - pip[i, j, k]) / dz
        )

    return gradient
end
