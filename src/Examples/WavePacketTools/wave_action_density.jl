# src/Examples/WavePacketTools/wave_action_density.jl

"""
```julia
wave_action_density(state::State, parameters::NamedTuple, x::Real, y::Real, z::Real)::Real
```

Return the wave-action density at ``\\left(x, y, z\\right)``.

The wave-action density is given by

```math
\\mathcal{A} = \\frac{\\bar{\\rho}}{2} \\frac{\\hat{\\omega} \\left|\\boldsymbol{k}\\right|^2}{N^4 \\left(k^2 + l^2\\right)} \\left|b_\\mathrm{w}\\right|^2,
```

where ``\\boldsymbol{k} = \\left(k, l, m\\right)^\\mathrm{T}`` is taken from `parameters` and ``\\bar{\\rho}``, ``\\hat{\\omega}``, ``N^2`` and ``b_\\mathrm{w}`` are calculated with `rhobar`, `omega`, `n2` and `bhat`, respectively.

# Arguments

  - `state`: Auxiliary model state.

  - `parameters`: Parameters of the wave-packet configuration.

  - `x`: Zonal position.

  - `y`: Meridional position.

  - `z`: Vertical position.

# See also

  - [`PinCFlow.Examples.WavePacketTools.bhat`](@ref)

  - [`PinCFlow.Examples.WavePacketTools.n2`](@ref)

  - [`PinCFlow.Examples.WavePacketTools.omega`](@ref)

  - [`PinCFlow.Examples.WavePacketTools.rhobar`](@ref)
"""
function wave_action_density end

function wave_action_density(
    state::State,
    parameters::NamedTuple,
    x::Real,
    y::Real,
    z::Real,
)::Real
    (; k, l, m) = parameters

    return n2(state, x, y, z) == 0.0 ? 0.0 :
           rhobar(state, x, y, z) / 2 *
           omega(state, parameters, x, y, z) *
           (k^2 + l^2 + m^2) / n2(state, x, y, z)^2 / (k^2 + l^2) *
           bhat(state, parameters, x, y, z)^2
end
