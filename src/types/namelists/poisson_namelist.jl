struct PoissonNamelist{A <: AbstractFloat, B <: Integer, C <: Bool}
  tolpoisson::A
  maxiterpoisson::B
  preconditioner::C
  dtau::A
  maxiteradi::B
  initialcleaning::C
  relative_tolerance::C
end

function PoissonNamelist(;
  tolpoisson = 1.0E-8,
  maxiterpoisson = 1000,
  preconditioner = true,
  dtau = 1.0E+0,
  maxiteradi = 2,
  initialcleaning = true,
  relative_tolerance = false,
)
  return PoissonNamelist(
    tolpoisson,
    maxiterpoisson,
    preconditioner,
    dtau,
    maxiteradi,
    initialcleaning,
    relative_tolerance,
  )
end
