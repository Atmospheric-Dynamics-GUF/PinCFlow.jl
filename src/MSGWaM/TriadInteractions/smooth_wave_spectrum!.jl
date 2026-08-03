function smooth_wave_spectrum! end

function smooth_wave_spectrum!(state::State)
    (; x_size, y_size) = state.namelists.domain
    (; filter_type, filter_order) = state.namelists.wkb

    (; spec_tend) = state.wkb
    (; wavespectrum, was_pred) = spec_tend
    (; kp, m) = spec_tend.spec_grid

    (; nxx, nyy, nzz) = state.domain
    (; jac) = state.grid

    if !(1 <= filter_order <= 4)
        error(
            "The wave-spectrum Shapiro-filter order must be between 1 and 4. ",
            "Received filter_order = ",
            filter_order,
        )
    end

    if !(filter_type isa Shapiro)
        error(
            "Wave-spectrum smoothing is currently implemented only ",
            "for the Shapiro filter.",
        )
    end

    # Convert N to the Jacobian-weighted density Q = J*N.
    @ivy for mi in eachindex(m), kpi in eachindex(kp),
        k in 1:nzz, j in 1:nyy, i in 1:nxx

        wavespectrum[i, j, k, kpi, mi] *= jac[i, j, k]
    end

    if x_size == y_size == 1
        apply_wave_spectrum_filter!(state, filter_type, Z())
    elseif x_size == 1
        apply_wave_spectrum_filter!(state, filter_type, YZ())
    elseif y_size == 1
        apply_wave_spectrum_filter!(state, filter_type, Z())  #for now the filter in only z direction have applied, later it should be replace by XZ()
    else
        apply_wave_spectrum_filter!(state, filter_type, XYZ())
    end

    # Convert the filtered J*N back to wave-action density N.
    @ivy for mi in eachindex(m), kpi in eachindex(kp),
        k in 1:nzz, j in 1:nyy, i in 1:nxx

        wavespectrum[i, j, k, kpi, mi] /= jac[i, j, k]
    end

    # Store the smoothed spectrum as the pre-interaction spectrum.
    copyto!(was_pred, wavespectrum)

    return nothing
end
