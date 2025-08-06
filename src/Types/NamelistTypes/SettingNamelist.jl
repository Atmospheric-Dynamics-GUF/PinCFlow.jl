"""
```julia
SettingNamelist{
    A <: AbstractModel,
    B <: AbstractTestCase,
    C <: AbstractBoundaries,
}
```

Namelist for parameters describing the general model setting.

```julia
SettingNamelist(;
    model::AbstractModel = PseudoIncompressible(),
    testcase::AbstractTestCase = MountainWave(),
    zboundaries::AbstractBoundaries = SolidWallBoundaries(),
)
```

Construct a `SettingNamelist` instance with the given keyword arguments as properties.

# Fields/Keywords

  - `model::A`: Dynamic equations.
  - `testcase::B`: Test case on wich the current simulation is based.
  - `zboundaries::C`: Vertical boundary conditions.
"""
struct SettingNamelist{
    A <: AbstractModel,
    B <: AbstractTestCase,
    C <: AbstractBoundaries,
}
    model::A
    testcase::B
    zboundaries::C
end

function SettingNamelist(;
    model::AbstractModel = PseudoIncompressible(),
    testcase::AbstractTestCase = MountainWave(),
    zboundaries::AbstractBoundaries = SolidWallBoundaries(),
)
    return SettingNamelist(model, testcase, zboundaries)
end
