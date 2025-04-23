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
