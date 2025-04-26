function set_boundary_rays!(state::State)
    (; testcase) = state.namelists.setting
    set_boundary_rays!(state, testcase)
    return
end

function set_boundary_rays!(state::State, testcase::AbstractTestCase)
    return
end

function set_boundary_rays!(state::State, testcase::AbstractWKBTestCase)
    (; wkb_mode) = state.namelists.wkb
    set_boundary_rays!(state, wkb_mode)
    return
end

function set_boundary_rays!(state::State, wkb_mode::SteadyState)
    (; sizex, sizey) = state.namelists.domain

    if sizex > 1
        set_zonal_boundary_rays!(state)
    end
    if sizey > 1
        set_meridional_boundary_rays!(state)
    end

    return
end

function set_boundary_rays!(
    state::State,
    wkb_mode::AbstractWKBMode,
)
    (; sizex, sizey) = state.namelists.domain
    (; zboundaries) = state.namelists.setting

    if sizex > 1
        set_zonal_boundary_rays!(state)
    end
    if sizey > 1
        set_meridional_boundary_rays!(state)
    end
    set_vertical_boundary_rays!(state, zboundaries)

    return
end
