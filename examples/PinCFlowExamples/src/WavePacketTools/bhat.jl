function bhat end

function bhat(
    state::State,
    parameters::NamedTuple,
    x::Real,
    y::Real,
    z::Real,
)::Real
    (; a0, m) = parameters

    return a0 * n2(state, x, y, z) / m * envelope(parameters, x, y, z)
end
