"""
```julia
GridNamelist
```

Namelist for parameters describing the grid.

```julia
GridNamelist(;
    vertical_grid_stretching::Function = zhat -> zhat,
    resolved_topography::Function = (x, y, t) -> 0.0,
    unresolved_topography::Function = (alpha, x, y, t) -> (0.0, 0.0, 0.0),
)::GridNamelist
```

Construct a `GridNamelist` instance with the given keyword arguments as properties, converting them to meet the type constraints.

The topography (both its resolved part and its unresolved part) may have a time dependence. At the beginning of each time step, it is updated accordingly. If that update resulted in a change, the vertical grid and atmospheric background fields are recomputed, the predictands are linearly interpolated to the new grid, and MS-GWaM's ray volumes are assigned to the new cells.

# Fields/Keywords

  - `vertical_grid_stretching::FunctionWrapper{Float64, Tuple{Float64}}`: Function that defines the vertical grid stretching ``\\tilde{z} \\left(\\hat{z}\\right)``.

  - `resolved_topography::FunctionWrapper{Float64, NTuple{2, Float64}}`: Function that returns the resolved topography at a specified horizontal position and time.

  - `unresolved_topography::FunctionWrapper{NTuple{3, Float64}, Tuple{Int, Float64, Float64}}`: Function that returns a specified spectral mode of the unresolved topography at a specified horizontal position and time.
"""
struct GridNamelist
    vertical_grid_stretching::FunctionWrapper{Float64, Tuple{Float64}}
    resolved_topography::FunctionWrapper{Float64, NTuple{3, Float64}}
    unresolved_topography::FunctionWrapper{
        NTuple{3, Float64},
        Tuple{Int, Float64, Float64, Float64},
    }
end

function GridNamelist(;
    vertical_grid_stretching::Function = zhat -> zhat,
    resolved_topography::Function = (x, y, t) -> 0.0,
    unresolved_topography::Function = (alpha, x, y, t) -> (0.0, 0.0, 0.0),
)::GridNamelist
    return GridNamelist(
        vertical_grid_stretching,
        resolved_topography,
        unresolved_topography,
    )
end
