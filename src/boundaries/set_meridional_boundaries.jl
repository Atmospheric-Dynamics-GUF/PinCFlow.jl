function set_meridional_boundaries!(
  state::State,
  variables::BoundaryPredictands,
  boundaries::PeriodicBoundaries,
)
  (; namelists, domain) = state
  (; rho, rhop, u, v, w, pip) = state.variables.predictands

  set_meridional_boundaries_of_field!(rho, namelists, domain)
  set_meridional_boundaries_of_field!(rhop, namelists, domain)
  set_meridional_boundaries_of_field!(u, namelists, domain)
  set_meridional_boundaries_of_field!(v, namelists, domain)
  set_meridional_boundaries_of_field!(w, namelists, domain)
  set_meridional_boundaries_of_field!(pip, namelists, domain)

  return
end

function set_meridional_boundaries!(
  state::State,
  variables::BoundaryReconstructions,
  boundaries::PeriodicBoundaries,
)
  (; namelists, domain) = state
  (rhotilde, rhoptilde, utilde, vtilde, wtilde) =
    state.variables.reconstructions

  set_meridional_boundaries_of_field!(rhotilde, namelists, domain)
  set_meridional_boundaries_of_field!(rhoptilde, namelists, domain)
  set_meridional_boundaries_of_field!(utilde, namelists, domain)
  set_meridional_boundaries_of_field!(vtilde, namelists, domain)
  set_meridional_boundaries_of_field!(wtilde, namelists, domain)

  return
end
