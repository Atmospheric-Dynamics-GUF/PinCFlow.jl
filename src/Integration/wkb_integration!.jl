"""
```julia
wkb_integration!(state::State, dtstage::AbstractFloat)
```
"""
function wkb_integration! end

function wkb_integration!(state::State, dtstage::AbstractFloat)
    (; nstages) = state.time

    apply_saturation_scheme!(state, dtstage)

    for rkstage in 1:nstages
        propagate_rays!(state, dtstage, rkstage)
    end

    split_rays!(state)
    shift_rays!(state)
    merge_rays!(state)
    set_boundary_rays!(state)

    compute_mean_flow_effect!(state)

    return
end
