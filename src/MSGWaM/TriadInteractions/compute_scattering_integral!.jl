function compute_scattering_integral! end

function compute_scattering_integral!(
    state::State,
    ii::Integer,
    jj::Integer,
    kk::Integer,
    triad_mode::Triad2D,
)
    (; x_size) = state.namelists.domain

    if x_size == 1
        compute_scattering_integral_discrete!(state, ii, jj, kk, triad_mode)
    else
        compute_scattering_integral_continuous!(state, ii, jj, kk, triad_mode)
    end

    return nothing
end