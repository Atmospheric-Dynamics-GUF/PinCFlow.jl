function set_boundary_ray_volumes!(state::State)
    (; wkb_mode) = state.namelists.wkb
    set_boundary_ray_volumes!(state, wkb_mode)
    return
end

function set_boundary_ray_volumes!(state::State, wkb_mode::SteadyState)
    return
end

function set_boundary_ray_volumes!(state::State, wkb_mode::SingleColumn)
    (; zboundaries) = state.namelists.setting
    set_vertical_boundary_ray_volumes!(state, zboundaries)
    return
end

function set_boundary_ray_volumes!(state::State, wkb_mode::MultiColumn)
    (; sizex, sizey) = state.namelists.domain
    (; zboundaries) = state.namelists.setting

    if sizex > 1
        set_zonal_boundary_ray_volumes!(state)
    end
    if sizey > 1
        set_meridional_boundary_ray_volumes!(state)
    end
    set_vertical_boundary_ray_volumes!(state, zboundaries)

    return
end
