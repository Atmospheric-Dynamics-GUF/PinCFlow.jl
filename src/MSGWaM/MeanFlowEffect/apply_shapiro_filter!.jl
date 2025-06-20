"""
    apply_shapiro_filter!(output::AbstractVector{<:AbstractFloat}, input::AbstractVector{<:AbstractFloat}, bounds::NTuple{2, <:Integer}, order::Val{1})

# Arguments

  - `output::AbstractVector{<:AbstractFloat}`: Output vector to store filtered values
  - `input::AbstractVector{<:AbstractFloat}`: Input vector to filter
  - `bounds::NTuple{2, <:Integer}`: Index bounds for filtering (start, end)
  - `order::Val{1}`: First-order Shapiro filter
"""
function apply_shapiro_filter!(
    output::AbstractVector{<:AbstractFloat},
    input::AbstractVector{<:AbstractFloat},
    bounds::NTuple{2, <:Integer},
    order::Val{1},
)
    for i in bounds[1]:bounds[2]
        output[i] = (input[i - 1] + input[i + 1] + 2 * input[i]) / 4
    end
    return
end

"""
    apply_shapiro_filter!(output::AbstractVector{<:AbstractFloat}, input::AbstractVector{<:AbstractFloat}, bounds::NTuple{2, <:Integer}, order::Val{2})

# Arguments

  - `output::AbstractVector{<:AbstractFloat}`: Output vector to store filtered values
  - `input::AbstractVector{<:AbstractFloat}`: Input vector to filter
  - `bounds::NTuple{2, <:Integer}`: Index bounds for filtering (start, end)
  - `order::Val{2}`: Second-order Shapiro filter
"""
function apply_shapiro_filter!(
    output::AbstractVector{<:AbstractFloat},
    input::AbstractVector{<:AbstractFloat},
    bounds::NTuple{2, <:Integer},
    order::Val{2},
)
    for i in bounds[1]:bounds[2]
        output[i] =
            (
                -input[i - 2] - input[i + 2] +
                4 * (input[i - 1] + input[i + 1]) +
                10 * input[i]
            ) / 16
    end
    return
end

"""
    apply_shapiro_filter!(output::AbstractVector{<:AbstractFloat}, input::AbstractVector{<:AbstractFloat}, bounds::NTuple{2, <:Integer}, order::Val{3})

# Arguments

  - `output::AbstractVector{<:AbstractFloat}`: Output vector to store filtered values
  - `input::AbstractVector{<:AbstractFloat}`: Input vector to filter
  - `bounds::NTuple{2, <:Integer}`: Index bounds for filtering (start, end)
  - `order::Val{3}`: Third-order Shapiro filter
"""
function apply_shapiro_filter!(
    output::AbstractVector{<:AbstractFloat},
    input::AbstractVector{<:AbstractFloat},
    bounds::NTuple{2, <:Integer},
    order::Val{3},
)
    for i in bounds[1]:bounds[2]
        output[i] =
            (
                input[i - 3] + input[i + 3] -
                6 * (input[i - 2] + input[i + 2]) +
                15 * (input[i - 1] + input[i + 1]) +
                44 * input[i]
            ) / 64
    end
    return
end

"""
    apply_shapiro_filter!(output::AbstractVector{<:AbstractFloat}, input::AbstractVector{<:AbstractFloat}, bounds::NTuple{2, <:Integer}, order::Val{4})

# Arguments

  - `output::AbstractVector{<:AbstractFloat}`: Output vector to store filtered values
  - `input::AbstractVector{<:AbstractFloat}`: Input vector to filter
  - `bounds::NTuple{2, <:Integer}`: Index bounds for filtering (start, end)
  - `order::Val{4}`: Fourth-order Shapiro filter
"""
function apply_shapiro_filter!(
    output::AbstractVector{<:AbstractFloat},
    input::AbstractVector{<:AbstractFloat},
    bounds::NTuple{2, <:Integer},
    order::Val{4},
)
    for i in bounds[1]:bounds[2]
        output[i] =
            (
                -input[i - 4] - input[i + 4] +
                8 * (input[i - 3] + input[i + 3]) -
                28 * (input[i - 2] + input[i + 2]) +
                56 * (input[i - 1] + input[i + 1]) +
                186 * input[i]
            ) / 256
    end
    return
end
