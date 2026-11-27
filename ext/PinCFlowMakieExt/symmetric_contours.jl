"""
```julia
symmetric_contours(
    minimum::Real,
    maximum::Real;
    colormap_name::Symbol = :seismic,
    digits::Integer = 2,
    number::Integer = 10,
)::Tuple{<:AbstractVector{<:Real}, <:Any}
```

Compute symmetric contour levels and return them and a correspondingly indexed colormap.

# Arguments

  - `minimum`: Smallest value to be plotted.

  - `maximum`: Largest value to be plotted.

# Keywords

    - `colormap_name`: Name under which the chosen colormap is registered.

    - `digits`: Number of digits of each level.

    - `number`: Number of contour levels.
"""
symmetric_contours

@ivy function symmetric_contours(
    minimum::Real,
    maximum::Real;
    colormap_name::Symbol = :seismic,
    digits::Integer = 2,
    number::Integer = 10,
)::Tuple{<:AbstractVector{<:Real}, <:Any}

    # Get the colormap and adjust the number of levels if necessary.
    colormap = to_colormap(colormap_name)
    (number - 1) % 2 != length(colormap) % 2 && (number += 1)

    # Add a tiny bit of noise if necessary.
    if minimum == maximum
        minimum -= eps(minimum)
        maximum += eps(maximum)
    elseif minimum == 0
        minimum -= eps(minimum)
    elseif maximum == 0
        maximum += eps(maximum)
    end

    # Determine the center.
    if sign(minimum) == sign(maximum)
        center = round((minimum + maximum) / 2; digits)
    else
        center = 0.0
    end

    # Determine the number of levels to the left and to the right.
    left_number =
        round(
            (center - minimum) / (maximum - minimum) * (number - 1) -
            ((number - 1) % 2) / 2,
        ) + ((number - 1) % 2) / 2
    right_number = number - 1 - left_number

    # Calculate the increment.
    increment = max(
        ceil((center - minimum) / left_number; digits),
        ceil((maximum - center) / right_number; digits),
    )
    (number - 1) % 2 != 0 &&
        increment / 2 != round(increment / 2; digits) &&
        (increment += 10.0^(-digits))

    # Set the levels.
    minimum = center - left_number * increment
    maximum = center + right_number * increment
    levels = round.(LinRange(minimum, maximum, number); digits)

    # Determine the indices for the colormap.
    midpoints = levels[1:(number - 1)] .+ (levels[2] .- levels[1]) ./ 2
    peak = max(abs(levels[1]), abs(levels[end]))
    indices = ceil.(Int, (1 .+ midpoints ./ peak) .* length(colormap) ./ 2)

    (levels, colormap[indices])
end
