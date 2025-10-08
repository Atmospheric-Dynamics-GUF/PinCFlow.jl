"""
```julia
SettingNamelist{A <: AbstractModel, B <: AbstractTestCase}
```

Namelist for parameters describing the general model setting.

```julia
SettingNamelist(;
    model::AbstractModel = PseudoIncompressible(),
    test_case::AbstractTestCase = MountainWave(),
)::SettingNamelist
```

Construct a `SettingNamelist` instance with the given keyword arguments as properties.

# Fields/Keywords

  - `model::A`: Dynamic equations.

  - `test_case::B`: Test case on which the current simulation is based.
"""
struct SettingNamelist{A <: AbstractModel, B <: AbstractTestCase}
    model::A
    test_case::B
end

function SettingNamelist(;
    model::AbstractModel = PseudoIncompressible(),
    test_case::AbstractTestCase = MountainWave(),
)::SettingNamelist
    return SettingNamelist(model, test_case)
end
