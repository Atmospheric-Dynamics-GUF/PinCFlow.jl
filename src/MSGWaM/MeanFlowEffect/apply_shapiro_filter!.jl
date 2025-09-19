"""
```julia
apply_shapiro_filter!(
    output::AbstractVector{<:AbstractFloat},
    input::AbstractVector{<:AbstractFloat},
    scope::UnitRange{<:Integer},
    order::Val{1},
)
```

Apply the first-order Shapiro filter to `input`.

The elements of `output` are given by

```math
\\widetilde{\\phi}_i = \\frac{1}{4} \\left(\\phi_{i - 1} + 2 \\phi_i + \\phi_{i + 1}\\right),
```

where ``\\phi_i`` are the elements of `input`.

```julia
apply_shapiro_filter!(
    output::AbstractVector{<:AbstractFloat},
    input::AbstractVector{<:AbstractFloat},
    scope::UnitRange{<:Integer},
    order::Val{2},
)
```

Apply the second-order Shapiro filter to `input`.

The elements of `output` are given by

```math
\\widetilde{\\phi}_i = \\frac{1}{16} \\left(- \\phi_{i - 2} + 4 \\phi_{i - 1} + 10 \\phi_i + 4 \\phi_{i + 1} - \\phi_{i + 2}\\right),
```

where ``\\phi_i`` are the elements of `input`.

```julia
apply_shapiro_filter!(
    output::AbstractVector{<:AbstractFloat},
    input::AbstractVector{<:AbstractFloat},
    scope::UnitRange{<:Integer},
    order::Val{3},
)
```

Apply the third-order Shapiro filter to `input`.

The elements of `output` are given by

```math
\\widetilde{\\phi}_i = \\frac{1}{64} \\left(\\phi_{i - 3} - 6 \\phi_{i - 2} + 15 \\phi_{i - 1} + 44 \\phi_i + 15 \\phi_{i + 1} - 6 \\phi_{i + 2} + \\phi_{i + 3}\\right),
```

where ``\\phi_i`` are the elements of `input`.

```julia
apply_shapiro_filter!(
    output::AbstractVector{<:AbstractFloat},
    input::AbstractVector{<:AbstractFloat},
    scope::UnitRange{<:Integer},
    order::Val{4},
)
```

Apply the fourth-order Shapiro filter to `input`.

The elements of `output` are given by

```math
\\widetilde{\\phi}_i = \\frac{1}{256} \\left(- \\phi_{i - 4} + 8 \\phi_{i - 3} - 28 \\phi_{i - 2} + 56 \\phi_{i - 1} + 186 \\phi_i + 56 \\phi_{i + 1} - 28 \\phi_{i + 2} + 8 \\phi_{i + 3} - \\phi_{i + 4}\\right),
```

where ``\\phi_i`` are the elements of `input`.

# Arguments

  - `output`: Filtered output vector.

  - `input`: Input vector.

  - `scope`: Index range.

  - `order`: Order of the Shapiro filter.
"""
function apply_shapiro_filter! end

function apply_shapiro_filter!(
    output::AbstractVector{<:AbstractFloat},
    input::AbstractVector{<:AbstractFloat},
    scope::UnitRange{<:Integer},
    order::Val{1},
)
    @ivy for i in scope
        output[i] = (input[i - 1] + input[i + 1] + 2 * input[i]) / 4
    end
    return
end

function apply_shapiro_filter!(
    output::AbstractVector{<:AbstractFloat},
    input::AbstractVector{<:AbstractFloat},
    scope::UnitRange{<:Integer},
    order::Val{2},
)
    @ivy for i in scope
        output[i] =
            (
                -input[i - 2] - input[i + 2] +
                4 * (input[i - 1] + input[i + 1]) +
                10 * input[i]
            ) / 16
    end
    return
end

function apply_shapiro_filter!(
    output::AbstractVector{<:AbstractFloat},
    input::AbstractVector{<:AbstractFloat},
    scope::UnitRange{<:Integer},
    order::Val{3},
)
    @ivy for i in scope
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

function apply_shapiro_filter!(
    output::AbstractVector{<:AbstractFloat},
    input::AbstractVector{<:AbstractFloat},
    scope::UnitRange{<:Integer},
    order::Val{4},
)
    @ivy for i in scope
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
