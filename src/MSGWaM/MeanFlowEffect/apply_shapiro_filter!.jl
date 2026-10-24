"""
```julia
apply_shapiro_filter!(
    output::AbstractVector{<:AbstractFloat},
    input::AbstractVector{<:AbstractFloat},
    order::Val{1},
    i::Integer,
)::Nothing
```

Apply the first-order Shapiro filter to `input` at the index `i`.

The elements of `output` are given by

```math
\\tilde{\\phi}_i = \\frac{1}{4} \\left(\\phi_{i - 1} + 2 \\phi_i + \\phi_{i + 1}\\right),
```

where ``\\phi_i`` are the elements of `input`.

```julia
apply_shapiro_filter!(
    output::AbstractVector{<:AbstractFloat},
    input::AbstractVector{<:AbstractFloat},
    order::Val{2},
    i::Integer,
)::Nothing
```

Apply the second-order Shapiro filter to `input` at the index `i`.

The elements of `output` are given by

```math
\\tilde{\\phi}_i = \\frac{1}{16} \\left(- \\phi_{i - 2} + 4 \\phi_{i - 1} + 10 \\phi_i + 4 \\phi_{i + 1} - \\phi_{i + 2}\\right),
```

where ``\\phi_i`` are the elements of `input`.

```julia
apply_shapiro_filter!(
    output::AbstractVector{<:AbstractFloat},
    input::AbstractVector{<:AbstractFloat},
    order::Val{3},
    i::Integer,
)::Nothing
```

Apply the third-order Shapiro filter to `input` at the index `i`.

The elements of `output` are given by

```math
\\tilde{\\phi}_i = \\frac{1}{64} \\left(\\phi_{i - 3} - 6 \\phi_{i - 2} + 15 \\phi_{i - 1} + 44 \\phi_i + 15 \\phi_{i + 1} - 6 \\phi_{i + 2} + \\phi_{i + 3}\\right),
```

where ``\\phi_i`` are the elements of `input`.

```julia
apply_shapiro_filter!(
    output::AbstractVector{<:AbstractFloat},
    input::AbstractVector{<:AbstractFloat},
    order::Val{4},
    i::Integer,
)::Nothing
```

Apply the fourth-order Shapiro filter to `input` at the index `i`.

The elements of `output` are given by

```math
\\tilde{\\phi}_i = \\frac{1}{256} \\left(- \\phi_{i - 4} + 8 \\phi_{i - 3} - 28 \\phi_{i - 2} + 56 \\phi_{i - 1} + 186 \\phi_i + 56 \\phi_{i + 1} - 28 \\phi_{i + 2} + 8 \\phi_{i + 3} - \\phi_{i + 4}\\right),
```

where ``\\phi_i`` are the elements of `input`.

# Arguments

  - `output`: Filtered output vector.

  - `input`: Input vector.

  - `order`: Order of the Shapiro filter.

  - `i`: Index.
"""
function apply_shapiro_filter! end

@ivy function apply_shapiro_filter!(
    output::AbstractVector{<:AbstractFloat},
    input::AbstractVector{<:AbstractFloat},
    order::Val{1},
    i::Integer,
)::Nothing
    output[i] = (input[i - 1] + input[i + 1] + 2 * input[i]) / 4
    nothing
end

@ivy function apply_shapiro_filter!(
    output::AbstractVector{<:AbstractFloat},
    input::AbstractVector{<:AbstractFloat},
    order::Val{2},
    i::Integer,
)::Nothing
    output[i] =
        (
            -input[i - 2] - input[i + 2] +
            4 * (input[i - 1] + input[i + 1]) +
            10 * input[i]
        ) / 16
    nothing
end

@ivy function apply_shapiro_filter!(
    output::AbstractVector{<:AbstractFloat},
    input::AbstractVector{<:AbstractFloat},
    order::Val{3},
    i::Integer,
)::Nothing
    output[i] =
        (
            input[i - 3] + input[i + 3] - 6 * (input[i - 2] + input[i + 2]) +
            15 * (input[i - 1] + input[i + 1]) +
            44 * input[i]
        ) / 64
    nothing
end

@ivy function apply_shapiro_filter!(
    output::AbstractVector{<:AbstractFloat},
    input::AbstractVector{<:AbstractFloat},
    order::Val{4},
    i::Integer,
)::Nothing
    output[i] =
        (
            -input[i - 4] - input[i + 4] + 8 * (input[i - 3] + input[i + 3]) -
            28 * (input[i - 2] + input[i + 2]) +
            56 * (input[i - 1] + input[i + 1]) +
            186 * input[i]
        ) / 256
    nothing
end
