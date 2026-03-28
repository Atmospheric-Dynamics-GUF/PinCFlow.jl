"""
```julia
apply_elastic_mode_selection!(
    state::State,
    field::AbstractVector{<:AbstractFloat},
)::Tuple{<:Integer, <:AbstractFloat}
```

Apply elastic mode selection to `field` by setting entries that are to be discarded to zero.

This method first determines the indices ``j`` that would sort `field` in descending order of its absolute values, so that one has

```math
\\left|E_{j_k}\\right| \\geq \\left|E_{j_{k + 1}}\\right| \\quad \\forall k,
```

where ``E`` represents `field`. Then, it calculates the effective participation ratio

```math
N = k_{\\max}^{- 1} \\min \\left[k_{\\max}, \\frac{\\left(\\sum_k \\left|E_{j_k}\\right|\\right)^2}{\\sum_k E_{j_k}^2}\\right],
```

where ``k_{\\max}`` is the minimum between the namelist parameter `maximum_mode_count` and the number of nonzero entries in `field`. Next, the neighbor similarity

```math
S = \\left(k_{\\max} - 1\\right)^{- 1} \\sum\\limits_{k = 1}^{k_{\\max} - 1} \\exp \\left(- \\frac{\\left|E_{j_k} / E_{j_{k + 1}}\\right| - 1}{\\delta}\\right),
```

is determined, where ``\\delta = 0.2`` defines one of three internal parameters. These two measures are combined in

```math
C = w_N N + w_S S,
```

where ``w_N = w_S = 0.5`` defines the other two internal parameters. Finally, the elastic mode count is calculated, following

```math
k_{\\star} = \\min \\left\\{k_{\\max}, \\max \\left[k_{\\min}, \\min \\left(k: F_k \\geq \\alpha\\right)\\right]\\right\\},
```

where

```math
F_k = \\frac{\\sum_{l = 1}^k \\left|E_{j_l}\\right|}{\\sum_l \\left|E_{j_l}\\right|}
```

defines the retained power fraction and

```math
\\alpha = \\alpha_{\\min} + \\left(\\alpha_{\\max} - \\alpha_{\\min}\\right) C
```

the target power fraction, with ``\\alpha_{\\min}`` and ``\\alpha_{\\max}`` representing the namelist parameters `minimum_power_fraction` and `maximum_power_fraction`, respectively. The elements at ``j_{k > k_{\\star}}`` are then set to zero and ``k_{\\star}`` and ``F_{k_{\\star}}`` are returned.

# Arguments

  - `state`: Model state.

  - `field`: Field to be modified.
"""
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
