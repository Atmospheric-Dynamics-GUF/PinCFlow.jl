function cartesian_indices end

function cartesian_indices(input::Tuple)::CartesianIndices
    return CartesianIndices(input)
end

function cartesian_indices(input::Tuple{CartesianIndices})::CartesianIndices
    return input[1]
end
