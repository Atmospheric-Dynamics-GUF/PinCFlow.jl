abstract type AbstractBackground end
struct Isothermal <: AbstractBackground end

abstract type AbstractCoriolis end
struct ConstantCoriolis <: AbstractCoriolis end