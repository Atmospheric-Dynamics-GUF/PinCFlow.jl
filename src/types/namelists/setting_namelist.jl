abstract type AbstractModel end
struct PseudoIncompressible <: AbstractModel end

abstract type AbstractTestCase end
struct MountainWave <: AbstractTestCase end

abstract type AbstractBoundaries end
struct PeriodicBoundaries <: AbstractBoundaries end
struct SolidWallBoundaries <: AbstractBoundaries end

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
  model = PseudoIncompressible(),
  testcase = MountainWave(),
  zboundaries = SolidWallBoundaries(),
)
  return SettingNamelist(model, testcase, zboundaries)
end
