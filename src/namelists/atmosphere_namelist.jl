
abstract type AbstractBackground end
struct Isothermal <: AbstractBackground end

abstract type AbstractCoriolis end
struct ConstantCoriolis <: AbstractCoriolis end

struct AtmosphereNamelist{
    A<:Bool,
    B<:AbstractFloat,
    C<:AbstractBackground,
    D<:Vector{<:AbstractFloat},
    E<:AbstractCoriolis,
}
    specifyreynolds::A
    reinv::B
    mu_viscous_dim::B
    background::C
    temp0_dim::B
    press0_dim::B
    backgroundflow_dim::D
    f_coriolis_dim::B
    corset::E
end

function AtmosphereNamelist(;
  specifyreynolds = false,
  reinv = 0.0,
  mu_viscous_dim = 0.0,
  background = Isothermal(),
  temp0_dim = 300.0,
  press0_dim = 100000.0,
  backgroundflow_dim = [0.0, 0.0, 0.0],
  f_coriolis_dim = 0.0,
  corset = ConstantCoriolis(),
)
  return AtmosphereNamelist(
    specifyreynolds,
    reinv,
    mu_viscous_dim,
    background,
    temp0_dim,
    press0_dim,
    backgroundflow_dim,
    f_coriolis_dim,
    corset,
  )
end
