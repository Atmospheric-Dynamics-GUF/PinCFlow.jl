function set_zonal_boundaries!(state::State, variables::BoundaryPredictands)
  (; namelists, domain) = state
  (; rho, rhop, u, v, w, pip) = state.variables.predictands

  set_all_zonal_boundary_layers!(rho, namelists, domain)
  set_all_zonal_boundary_layers!(rhop, namelists, domain)
  set_all_zonal_boundary_layers!(u, namelists, domain)
  set_all_zonal_boundary_layers!(v, namelists, domain)
  set_all_zonal_boundary_layers!(w, namelists, domain)
  set_all_zonal_boundary_layers!(pip, namelists, domain)

  return
end

function set_zonal_boundaries!(state::State, variables::BoundaryReconstructions)
  (; namelists, domain) = state
  (; rhotilde, rhoptilde, utilde, vtilde, wtilde) =
    state.variables.reconstructions

  set_all_zonal_boundary_layers!(rhotilde, namelists, domain)
  set_all_zonal_boundary_layers!(rhoptilde, namelists, domain)
  set_all_zonal_boundary_layers!(utilde, namelists, domain)
  set_all_zonal_boundary_layers!(vtilde, namelists, domain)
  set_all_zonal_boundary_layers!(wtilde, namelists, domain)

  return
end
