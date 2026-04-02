function thetabar end

function thetabar(state::State, x::Real, y::Real, z::Real)::Real
    (; atmosphere) = state
    (; thetaref) = state.constants

    return atmosphere.thetabar[ijk(state, x, y, z)] .* thetaref
end
