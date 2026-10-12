# src/Examples/WavePacketTools/omega.jl

"""
```julia
omega(
    state::State,
    parameters::NamedTuple,
    x::Real,
    y::Real,
    z::Real,
)::Real
```

Return the intrinsic frequency at ``\\left(x, y, z\\right)``.

Consistent with the dispersion relation, the intrinsic frequency is given by

```math
\\hat{\\omega} = - \\sqrt{\\frac{N^2 \\left(k^2 + l^2\\right) + f^2 m^2}{\\left|\\boldsymbol{k}\\right|^2}},
```

where ``\\boldsymbol{k} = \\left(k, l, m\\right)^\\mathrm{T}`` and ``f`` are taken from `parameters` and ``N^2`` is calculated with `n2`.

# Arguments

  - `state`: Auxiliary model state (needed by `n2`).

  - `parameters`: Parameters of the wave-packet configuration.

  - `x`: Zonal position.

  - `y`: Meridional position.

  - `z`: Vertical position.

# See also

  - [`PinCFlow.Examples.WavePacketTools.n2`](@ref)
"""
function omega end

function omega(
    state::State,
    parameters::NamedTuple,
    x::Real,
    y::Real,
    z::Real,
)::Real
    (; coriolis_frequency) = state.namelists.atmosphere
    (; k, l, m) = parameters

    return -sqrt(
        (n2(state, x, y, z) * (k^2 + l^2) + coriolis_frequency^2 * m^2) /
        (k^2 + l^2 + m^2),
    )
end
