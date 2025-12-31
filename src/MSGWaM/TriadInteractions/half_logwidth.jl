"""
```julia
half_logwidth(
    state::State,
    xr::AbstractFloat,
    yr::AbstractFloat,
    dxr::AbstractFloat,
    dyr::AbstractFloat,
)::NTuple{4, <:Integer}
```

From the given horizontal ray-volume position and extent, determine the indices of the grid cells that contain the ray-volume edges and return them (in the order left, right, backward and forward).

# Arguments

  - `state`: Model state.

  - `xr`: Ray-volume position in ``x``.

  - `yr`: Ray-volume position in ``y``.

  - `dxr`: Ray-volume extent in ``x``.

  - `dyr`: Ray-volume extent in ``y``.
"""
function half_logwidth end

function half_logwidth(
    logarray::AbstractVector{<:AbstractFloat},
)::NTuple{2, <:AbstractVector{<:AbstractFloat}}
    
    uper_half_width = zeros(length(logarray))
    lower_half_width = zeros(length(logarray))
    for i in eachindex(logarray)
        if i == 1
            uper_half_width[i] = (logarray[i+1] - logarray[i]) / 2
            lower_half_width[i] = (logarray[i+1] - logarray[i]) / 2
        elseif i == length(logarray)
            uper_half_width[i] = (logarray[i] - logarray[i-1]) / 2
            lower_half_width[i] = (logarray[i] - logarray[i-1]) / 2
        else
            uper_half_width[i] = (logarray[i+1] - logarray[i]) / 2
            lower_half_width[i] = (logarray[i] - logarray[i-1]) / 2
        end
    end
    
    return (uper_half_width, lower_half_width)
end
