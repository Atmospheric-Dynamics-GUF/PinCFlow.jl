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
