function n2 end

function n2(state::State, x::Real, y::Real, z::Real)::Real
    (; atmosphere) = state
    (; tref) = state.constants

    return atmosphere.n2[ijk(state, x, y, z)] ./ tref .^ 2
end
