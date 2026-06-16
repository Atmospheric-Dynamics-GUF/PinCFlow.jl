function compute_nearest_index end

function compute_nearest_index(k::AbstractVector{<:AbstractFloat}, 
    kpi::AbstractFloat)::Integer
    i = searchsortedfirst(k, kpi)

    if i <= 1
        return -1
    elseif i > length(k)
        return -1
    else
        return abs(k[i] - kpi) < abs(k[i-1] - kpi) ? i : i-1
    end
end