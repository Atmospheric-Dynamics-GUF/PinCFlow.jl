"""
```julia
@dispatch_input_file_format(input::Expr)
```

Macro that makes value dispatch static for the `input_file_format` parameter of `OutputNamelist`.

The parameter can take any of the following values:

  - `:HDF5`

  - `:NetCDF4`

# Arguments

  - `input`: Input expression with `Val` calls.
"""
macro dispatch_input_file_format end

macro dispatch_input_file_format(input::Expr)
    return esc(quote
        @dispatch (:HDF5, :NetCDF4) $(input)
    end)
end
