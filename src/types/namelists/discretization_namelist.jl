abstract type AbstractLimiter end
struct MCVariant <: AbstractLimiter end

struct DiscretizationNamelist{
  A <: AbstractFloat,
  B <: Bool,
  C <: AbstractLimiter,
}
  cfl::A
  dtmin_dim::A
  dtmax_dim::A
  adaptive_time_step::B
  limitertype::C
end

function DiscretizationNamelist(;
  cfl = 0.5,
  dtmin_dim = 1.0e-6,
  dtmax_dim = 1.0e3,
  adaptive_time_step = true,
  limitertype = MCVariant(),
)
  return DiscretizationNamelist(
    cfl,
    dtmin_dim,
    dtmax_dim,
    adaptive_time_step,
    limitertype,
  )
end
