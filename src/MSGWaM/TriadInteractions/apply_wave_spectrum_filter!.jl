function apply_wave_spectrum_filter! end

function apply_wave_spectrum_filter!(
    state::State,
    ::Shapiro,
    ::Z,
)
    (; spec_tend) = state
    (; wavespectrum, was_pred) = spec_tend
    (; kp, m) = spec_tend.spec_grid

    (; filter_order) = state.namelists.wkb
    (; nbz) = state.namelists.domain
    (; nxx, nyy, k0, k1) = state.domain

    if nbz < filter_order
        error(
            "Error in apply_wave_spectrum_filter!: ",
            "nbz < filter_order.",
        )
    end

    copyto!(was_pred, wavespectrum)

    @ivy for mi in eachindex(m), kpi in eachindex(kp),
        j in 1:nyy, i in 1:nxx

        apply_shapiro_filter!(
            wavespectrum[i, j, :, kpi, mi],
            was_pred[i, j, :, kpi, mi],
            k0:k1,
            Val(filter_order),
        )
    end

    return nothing
end

function apply_wave_spectrum_filter!(
    state::State,
    ::Shapiro,
    ::X,
)
    (; spec_tend) = state
    (; wavespectrum, was_pred) = spec_tend
    (; kp, m) = spec_tend.spec_grid

    (; filter_order) = state.namelists.wkb
    (; nbx) = state.namelists.domain
    (; nyy, nzz, i0, i1) = state.domain

    if nbx < filter_order
        error(
            "Error in apply_wave_spectrum_filter!: ",
            "nbx < filter_order.",
        )
    end

    copyto!(was_pred, wavespectrum)

    @ivy for mi in eachindex(m), kpi in eachindex(kp),
        k in 1:nzz, j in 1:nyy

        apply_shapiro_filter!(
            wavespectrum[:, j, k, kpi, mi],
            was_pred[:, j, k, kpi, mi],
            i0:i1,
            Val(filter_order),
        )
    end

    return nothing
end

function apply_wave_spectrum_filter!(
    state::State,
    ::Shapiro,
    ::Y,
)
    (; spec_tend) = state
    (; wavespectrum, was_pred) = spec_tend
    (; kp, m) = spec_tend.spec_grid

    (; filter_order) = state.namelists.wkb
    (; nby) = state.namelists.domain
    (; nxx, nzz, j0, j1) = state.domain

    if nby < filter_order
        error(
            "Error in apply_wave_spectrum_filter!: ",
            "nby < filter_order.",
        )
    end

    copyto!(was_pred, wavespectrum)

    @ivy for mi in eachindex(m), kpi in eachindex(kp),
        k in 1:nzz, i in 1:nxx

        apply_shapiro_filter!(
            wavespectrum[i, :, k, kpi, mi],
            was_pred[i, :, k, kpi, mi],
            j0:j1,
            Val(filter_order),
        )
    end

    return nothing
end

function apply_wave_spectrum_filter!(
    state::State,
    filter_type::Shapiro,
    ::XZ,
)
    apply_wave_spectrum_filter!(state, filter_type, X())
    apply_wave_spectrum_filter!(state, filter_type, Z())

    return nothing
end

function apply_wave_spectrum_filter!(
    state::State,
    filter_type::Shapiro,
    ::YZ,
)
    apply_wave_spectrum_filter!(state, filter_type, Y())
    apply_wave_spectrum_filter!(state, filter_type, Z())

    return nothing
end

function apply_wave_spectrum_filter!(
    state::State,
    filter_type::Shapiro,
    ::XYZ,
)
    apply_wave_spectrum_filter!(state, filter_type, X())
    apply_wave_spectrum_filter!(state, filter_type, Y())
    apply_wave_spectrum_filter!(state, filter_type, Z())

    return nothing
end