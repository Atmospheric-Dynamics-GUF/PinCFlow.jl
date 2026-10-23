# src/Examples/WavePacketTools/pihat.jl

"""
```julia
pihat(state::State, parameters::NamedTuple, x::Real, y::Real, z::Real)::Real
```

Return the Exner-pressure wave amplitude at ``\\left(x, y, z\\right)``.

Consistent with the polarization relations, the Exner-pressure wave amplitude is given by

```math
\\pi_\\mathrm{w} = \\frac{\\kappa}{R \\bar{\\theta}} \\frac{i}{m} \\frac{\\hat{\\omega}^2 - N^2}{N^2} b_\\mathrm{w},
```

where ``\\kappa`` and ``R`` are taken from `state.constants` and ``m`` from `parameters`, whereas ``\\bar{\\theta}``, ``N^2``, ``\\hat{\\omega}`` and ``b_\\mathrm{w}`` are calculated with `thetabar`, `n2`, `omega` and `bhat`, respectively.

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

  - [`PinCFlow.Examples.WavePacketTools.thetabar`](@ref)
"""
function pihat end

function pihat(
    state::State,
    parameters::NamedTuple,
    x::Real,
    y::Real,
    z::Real,
)::Number
    (; kappa, rsp) = state.constants
    (; m) = parameters

    return n2(state, x, y, z) == 0.0 ? 0.0 :
           kappa / rsp / thetabar(state, x, y, z) * 1im / m *
           (omega(state, parameters, x, y, z)^2 - n2(state, x, y, z)) /
           n2(state, x, y, z) * bhat(state, parameters, x, y, z)
end
