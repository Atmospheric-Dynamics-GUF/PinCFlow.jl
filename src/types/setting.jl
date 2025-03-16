abstract type AbstractModel end
struct PseudoIncompressible <: AbstractModel end

abstract type AbstractTestCase end
struct MountainWave <: AbstractTestCase end