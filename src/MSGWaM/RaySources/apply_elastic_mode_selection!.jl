function apply_elastic_mode_selection! end

function apply_elastic_mode_selection!(
    state::State,
    field::AbstractVector{<:AbstractFloat},
)::Tuple{<:Integer, <:AbstractFloat}
    (;
        elastic_mode_selection,
        minimum_mode_count,
        maximum_mode_count,
        minimum_power_fraction,
        maximum_power_fraction,
    ) = state.namelists.wkb
    (; sorted_wave_mode_indices) = state.wkb.elastic_mode_selection

    # Define the internal parameters of the algorithm.
    delta = 0.2
    wn = 0.5
    ws = 0.5

    # Get the indices of the sorted field.
    sortperm!(sorted_wave_mode_indices, field; lt = >=, by = abs)

    # Exclude modes that are zero.
    stop = findfirst(j -> field[j] == 0, sorted_wave_mode_indices)
    @ivy if stop === nothing
        j = sorted_wave_mode_indices
    elseif stop > 1
        j = sorted_wave_mode_indices[1:(stop - 1)]
    else
        return (0, 1.0)
    end

    # Determine the index bounds.
    kmin = min(minimum_mode_count, length(j))
    kmax = min(maximum_mode_count, length(j))

    # Compute the participation ratio.
    n = min(sum(abs, field)^2 / sum(a -> a^2, field), kmax) / kmax

    # Compute the neighbor similarity.
    s = 0
    @ivy for k in 1:(kmax - 1)
        s += exp(-(abs(field[j[k]]) / abs(field[j[k + 1]]) - 1) / delta)
    end
    s /= kmax - 1

    # Compute the combined measure.
    c = wn * n + ws * s

    # Compute the target power fraction.
    alpha =
        minimum_power_fraction +
        (maximum_power_fraction - minimum_power_fraction) * c

    # Determine the optimal number of modes.
    k = 0
    f = 0.0
    @ivy while k < maximum_mode_count && (k < minimum_mode_count || f < alpha)
        k += 1
        f += abs(field[j[k]]) / sum(abs, field)
    end

    # Discard the remaining modes.
    @ivy field[j[(k + 1):end]] .= 0

    return (k, f)
end
