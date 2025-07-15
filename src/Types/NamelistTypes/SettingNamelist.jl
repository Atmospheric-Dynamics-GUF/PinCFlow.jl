"""
```julia
SettingNamelist{
    A <: AbstractModel,
    B <: AbstractTestCase,
    C <: AbstractBoundaries,
}
```

Namelist for level of compressibility, test case and vertical boundary conditions (see constructor for parameter descriptions).
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

"""
```julia
SettingNamelist(;
    model = PseudoIncompressible(),
    testcase = MountainWave(),
    zboundaries = SolidWallBoundaries(),
)
```

Core simulation configuration specifying equation set, test case, and vertical boundary conditions.

# Arguments

  - `model`: Governing equation set. Options: `PseudoIncompressible()`, `Boussinesq()`, `Compressible()`
  - `testcase`: Initial/forcing configuration. Options: `MountainWave()`, `WKBMountainWave()`
  - `zboundaries`: Vertical boundary treatment. Options: `SolidWallBoundaries()`, `PeriodicBoundaries()`

# Returns

  - `SettingNamelist`: Configuration container for simulation setup
"""
function SettingNamelist(;
    model = PseudoIncompressible(),
    testcase = MountainWave(),
    zboundaries = SolidWallBoundaries(),
)
    return SettingNamelist(model, testcase, zboundaries)
end
