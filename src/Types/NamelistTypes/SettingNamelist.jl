"""
```julia
SettingNamelist{A <: AbstractModel, B <: AbstractTestCase}
```

Namelist for parameters describing the general model setting.

```julia
SettingNamelist(;
    model::AbstractModel = PseudoIncompressible(),
    testcase::AbstractTestCase = MountainWave(),
)::SettingNamelist
```

Construct a `SettingNamelist` instance with the given keyword arguments as properties.

# Fields/Keywords

  - `model::A`: Dynamic equations.

  - `testcase::B`: Test case on wich the current simulation is based.
"""
struct SettingNamelist{A <: AbstractModel, B <: AbstractTestCase}
    model::A
    testcase::B
end

function SettingNamelist(;
    model::AbstractModel = PseudoIncompressible(),
    testcase::AbstractTestCase = MountainWave(),
)::SettingNamelist
    return SettingNamelist(model, testcase)
end
