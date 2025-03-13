
struct PoissonNamelist{
    A<:AbstractFloat,
    B<:Integer,
    C<:Bool,
}
    tolpoisson::A
    maxiterpoisson::B
    preconditioner::C
    dtau::A
    maxiteradi::B
    initialcleaning::C
    correctmomentum::C
    relative_tolerance::C
end

function PoissonNamelist(;
  tolpoisson = 1.0e-8,
  maxiterpoisson = 1000,
  preconditioner = true,
  dtau = 4.0e-4,
  maxiteradi = 2,
  initialcleaning = true,
  correctmomentum = true,
  relative_tolerance = false,
)
  return PoissonNamelist(
    tolpoisson,
    maxiterpoisson,
    preconditioner,
    dtau,
    maxiteradi,
    initialcleaning,
    correctmomentum,
    relative_tolerance,
  )
end
