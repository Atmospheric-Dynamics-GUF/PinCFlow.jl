function compute_edges end

function compute_edges(logarray::AbstractVector{<:Real})::AbstractVector{<:Real}
    ll = length(logarray)
    if ll == 1
        return zeros(1)
    end
    edges = similar(logarray, promote_type(eltype(logarray), Float64), ll + 1)

    # Lower and upper outer edges (extrapolated half-step at ends)
    edges[1]     = logarray[1]   - (logarray[2]   - logarray[1])   / 2
    edges[end]   = logarray[end] + (logarray[end] - logarray[end-1]) / 2

    # Interior edges are midpoints between centers
    @inbounds for i in 2:ll
        edges[i] = (logarray[i-1] + logarray[i]) / 2
    end
    return edges
end