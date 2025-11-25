"""
```julia
compute_norms()::NTuple{2, NamedTuple}
```

Return a `NamedTuple` of ``L_2`` norms and a `NamedTuple` of ``L_\\infty`` norms, computed from the datasets in the first HDF5 file found in the current directory, before deleting that file.
"""
function compute_norms end

function compute_norms()::NTuple{2, NamedTuple}
    file = [name for name in readdir(".") if endswith(name, ".h5")][1]

    (l2, linf) = h5open(file, "r") do data
        l2 = NamedTuple(
            Symbol(key) => norm(read(data[key]), 2) for
            key in keys(data) if typeof(data[key]) <: HDF5.Dataset
        )
        linf = NamedTuple(
            Symbol(key) => norm(read(data[key]), Inf) for
            key in keys(data) if typeof(data[key]) <: HDF5.Dataset
        )
        return (l2, linf)
    end

    rm(file)

    return (l2, linf)
end
