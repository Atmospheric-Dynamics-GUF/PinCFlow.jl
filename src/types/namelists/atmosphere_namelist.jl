struct AtmosphereNamelist{
  A <: Bool,
  B <: AbstractFloat,
  C <: AbstractBackground,
  D <: AbstractVector{<:AbstractFloat},
  E <: AbstractCoriolis,
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
  reinv = 0.0E+0,
  mu_viscous_dim = 0.0E+0,
  background = Isothermal(),
  temp0_dim = 3.0E+2,
  press0_dim = 1.0E+5,
  backgroundflow_dim = [0.0E+0, 0.0E+0, 0.0E+0],
  f_coriolis_dim = 0.0E+0,
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
