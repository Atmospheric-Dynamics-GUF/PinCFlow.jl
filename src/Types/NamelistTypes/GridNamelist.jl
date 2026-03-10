"""
```julia
GridNamelist{A <: Float64, B <: Function, C <: Function}
```

Namelist for parameters describing the grid.

```julia
GridNamelist(;
    vertical_grid_stretching::Function = zhat -> zhat,
    resolved_topography::Function = (x, y) -> 0.0,
    unresolved_topography::Function = (alpha, x, y) -> (0.0, 0.0, 0.0),
)::GridNamelist
```

Construct a `GridNamelist` instance with the given keyword arguments as properties, converting them to meet the type constraints.

# Fields/Keywords

  - `vertical_grid_stretching::A`: Function that defines the vertical grid stretching ``\\tilde{z} \\left(\\hat{z}\\right)``.

  - `resolved_topography::B`: Function that returns the resolved topography at a specified horizontal position.

  - `resolved_topography::C`: Function that returns a specified spectral mode of the unresolved topography at a specified horizontal position.
"""
struct GridNamelist{A <: Function, B <: Function, C <: Function}
    vertical_grid_stretching::A
    resolved_topography::B
    unresolved_topography::C
end

function GridNamelist(;
    vertical_grid_stretching::Function = zhat -> zhat,
    resolved_topography::Function = (x, y) -> 0.0,
    unresolved_topography::Function = (alpha, x, y) -> (0.0, 0.0, 0.0),
)::GridNamelist
    return GridNamelist(
        vertical_grid_stretching,
        resolved_topography,
        unresolved_topography,
    )
end
