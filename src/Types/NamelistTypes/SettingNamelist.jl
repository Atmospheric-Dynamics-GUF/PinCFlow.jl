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
    SettingNamelist(; 
        model = PseudoIncompressible(),
        testcase = MountainWave(),
        zboundaries = SolidWallBoundaries()
    )

Core simulation configuration specifying equation set, test case, and vertical boundary conditions.

# Arguments

  - `model = PseudoIncompressible()`: Governing equation set. Options: `PseudoIncompressible()`, `Boussinesq()`, `Compressible()`
  - `testcase = MountainWave()`: Initial/forcing configuration. Options: `MountainWave()`, `WKBMountainWave()`
  - `zboundaries = SolidWallBoundaries()`: Vertical boundary treatment. Options: `SolidWallBoundaries()`, `PeriodicBoundaries()`

# Usage

Determines which equations are solved, how initial conditions are set, and boundary condition implementation.
"""
function SettingNamelist(;
    model = PseudoIncompressible(),
    testcase = MountainWave(),
    zboundaries = SolidWallBoundaries(),
)
    return SettingNamelist(model, testcase, zboundaries)
end
