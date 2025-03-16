struct SettingNamelist{A<:AbstractModel,B<:AbstractTestCase}
    model::A
    testcase::B
end

function SettingNamelist(;
  model = PseudoIncompressible(),
  testcase = MountainWave(),
)
  return SettingNamelist(model, testcase)
end
