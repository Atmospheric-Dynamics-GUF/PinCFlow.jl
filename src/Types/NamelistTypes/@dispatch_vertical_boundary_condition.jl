"""
```julia
@dispatch_vertical_boundary_condition(input::Expr)
```

Macro that makes value dispatch static for the `vertical_boundary_condition` parameter of `DomainNamelist`.

The parameter can take any of the following values:

  - `:SolidWall`

  - `:Periodic`

# Arguments

  - `input`: Input expression with `Val` calls.
"""
macro dispatch_vertical_boundary_condition end

macro dispatch_vertical_boundary_condition(input::Expr)
    return esc(
        quote
            @dispatch (:SolidWall, :Periodic) $(
                input
            )
        end,
    )
end
