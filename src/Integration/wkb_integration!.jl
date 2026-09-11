"""
```julia
wkb_integration!(state::State, dt::AbstractFloat)
```

Use MS-GWaM to update the properties of the unresolved gravity-wave field and compute its impact on the resolved flow.

In the first step, MS-GWaM's saturation scheme is applied to account for the impact of wave breaking on the wave-action field. Next, the ray volumes are propagated via integration of the ray equations with a Runge-Kutta time step. Ray volumes that have grown larger than the cells of the model grid are then split before their array positions are shifted such that they are attributed to the correct grid cells. Afterwards, ray volumes are merged in cells where their count exceeds a specified threshold. Finally, the ray volumes in boundary and halo cells are updated and the mean-flow impact is calculated.

# Arguments

  - `state`: Model state.

  - `dt`: Time step.
"""
function wkb_integration! end

function wkb_integration!(state::State, dt::AbstractFloat)
    (; nstages) = state.time

    apply_saturation_scheme!(state, dt)

    for rkstage in 1:nstages
        propagate_rays!(state, dt, rkstage)
    end

    split_rays!(state)
    shift_rays!(state)
    merge_rays!(state)
    set_boundary_rays!(state)

    compute_mean_flow_effect!(state, dt)

    return
end
