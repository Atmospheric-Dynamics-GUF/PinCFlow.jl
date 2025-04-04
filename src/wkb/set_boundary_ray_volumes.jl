function set_boundary_ray_volumes!(state::State)
  (; sizex, sizey) = state.namelists.domain
  (; zboundaries) = state.namelists.setting
  (; steady_state) = state.namelists.wkb

  if steady_state
    return
  else
    if sizex > 1
      set_zonal_boundary_ray_volumes!(state)
    end
    if sizey > 1
      set_meridional_boundary_ray_volumes!(state)
    end
    set_vertical_boundary_ray_volumes!(state, zboundaries)
    return
  end
end
