abstract type AbstractModel end
struct PseudoIncompressible <: AbstractModel end

abstract type AbstractTestCase end
struct MountainWave <: AbstractTestCase end

struct SettingNamelist{A <: AbstractModel, B <: AbstractTestCase}
  model::A
  testcase::B
end

function SettingNamelist(;
  model = PseudoIncompressible(),
  testcase = MountainWave(),
)
  return SettingNamelist(model, testcase)
end
