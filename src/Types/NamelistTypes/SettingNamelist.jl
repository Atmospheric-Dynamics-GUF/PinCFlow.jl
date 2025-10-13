"""
```julia
SettingNamelist{A <: AbstractModel}
```

Namelist for parameters describing the general model setting.

```julia
SettingNamelist(;
    model::AbstractModel = PseudoIncompressible(),
)::SettingNamelist
```

Construct a `SettingNamelist` instance with the given keyword arguments as properties.

# Fields/Keywords

  - `model::A`: Dynamic equations.
"""
struct SettingNamelist{A <: AbstractModel}
    model::A
end

function SettingNamelist(;
    model::AbstractModel = PseudoIncompressible(),
)::SettingNamelist
    return SettingNamelist(model)
end
