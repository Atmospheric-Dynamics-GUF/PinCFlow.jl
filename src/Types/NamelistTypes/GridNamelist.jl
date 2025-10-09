"""
```julia
GridNamelist{A <: AbstractFloat, B <: Function, C <: Function}
```

Namelist for parameters describing the grid.

```julia
GridNamelist(;
    stretch_exponent::AbstractFloat = 1.0E+0,
    resolved_topography::Function = (x, y) -> 0.0,
    unresolved_topography::Function = (alpha, x, y) -> (0.0, 0.0, 0.0),
)::GridNamelist
```

Construct a `GridNamelist` instance with the given keyword arguments as properties.

# Fields/Keywords

  - `stretch_exponent::A`: Vertical-grid-stretching parameter.

  - `resolved_topography::B`: Function that returns the resolved topography at a specified horizontal position.

  - `resolved_topography::C`: Function that returns a specified spectral mode of the unresolved topography at a specified horizontal position.
"""
struct GridNamelist{A <: AbstractFloat, B <: Function, C <: Function}
    stretch_exponent::A
    resolved_topography::B
    unresolved_topography::C
end

function GridNamelist(;
    stretch_exponent::AbstractFloat = 1.0E+0,
    resolved_topography::Function = (x, y) -> 0.0,
    unresolved_topography::Function = (alpha, x, y) -> (0.0, 0.0, 0.0),
)::GridNamelist
    return GridNamelist(
        stretch_exponent,
        resolved_topography,
        unresolved_topography,
    )
end
