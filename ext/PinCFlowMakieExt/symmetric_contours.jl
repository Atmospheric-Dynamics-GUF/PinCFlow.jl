function symmetric_contours(
    minimum::AbstractFloat,
    maximum::AbstractFloat;
    number::Integer = 10,
    colormap_name::Symbol = :seismic,
)::Tuple{<:LinRange{<:AbstractFloat, <:Integer}, <:Any}

    # Get the colormap and adjust the number of levels if necessary.
    colormap = to_colormap(colormap_name)
    (number - 1) % 2 != length(colormap) % 2 && (number += 1)

    # Compute contour levels.
    @ivy if minimum == -maximum ||
            sign(minimum) == sign(maximum) ||
            minimum == 0 ||
            maximum == 0
        levels = LinRange(minimum, maximum, number)
    else
        peak = max(abs(minimum), abs(maximum))
        factor = ceil(Int, 2 * peak / (maximum - minimum))
        if number % 2 > 0 && factor % 2 == 0
            factor += 1
        end
        levels = LinRange(-peak, peak, factor * number)
        if peak > -minimum
            while levels[end - number + 1] > minimum
                levels = LinRange(-peak, peak, length(levels) - 2)
            end
            levels = levels[(end - number + 1):end]
        elseif peak > maximum
            while levels[number] < maximum
                levels = LinRange(-peak, peak, length(levels) - 2)
            end
            levels = levels[1:number]
        end
    end

    # Determine the indices for the colormap.
    @ivy midpoints = levels[1:(number - 1)] .+ (levels[2] .- levels[1]) ./ 2
    @ivy peak = max(abs(levels[1]), abs(levels[end]))
    indices = ceil.(Int, (1 .+ midpoints ./ peak) .* length(colormap) ./ 2)

    return (levels, colormap[indices])
end
