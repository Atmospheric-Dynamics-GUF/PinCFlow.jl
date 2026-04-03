# examples/PinCFlowExamples/src/WavePacketTools/ijk.jl

function ijk(state::State, x::Real, y::Real, z::Real)::CartesianIndex
    (; lref) = state.constants
    (; grid) = state

    i = argmin(abs.(x .- grid.x .* lref))
    j = argmin(abs.(y .- grid.y .* lref))
    k = argmin(abs.(z .- grid.zc[i, j, :] .* lref))

    return CartesianIndex(i, j, k)
end
