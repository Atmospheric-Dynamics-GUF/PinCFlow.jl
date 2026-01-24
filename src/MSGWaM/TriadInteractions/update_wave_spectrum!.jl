function update_wave_spectrum! end

function update_wave_spectrum!(
    state::State,
    ii::Integer,
    jj::Integer,
    kk::Integer,
    dtau::AbstractFloat, 
    triad_mode::Union{Triad2D, Triad3DIso}
)
    (; spec_tend) = state.wkb
    (; wavespectrum, col_int) = spec_tend
    (; kp, m) = spec_tend.spec_grid

    compute_scattering_integral!(state, ii, jj, kk, triad_mode)


    #Euler method 
    """
    for mi in eachindex(m),
        kpi in eachindex(kp)

        if  col_int[kpi, mi] != 0
            wavespectrum[ii, jj, kk, kpi, mi] += dtau * col_int[kpi, mi]  ##Euler method
        end

    end
    """


    #RK2 method
    was_copy = spec_tend.wavespectrum
    for mi in eachindex(m),
        kpi in eachindex(kp)

        if  col_int[kpi, mi] != 0
            wavespectrum[ii, jj, kk, kpi, mi] += 0.5 * dtau * col_int[kpi, mi]  
        end

    end

    compute_scattering_integral!(state, ii, jj, kk, triad_mode)

    for mi in eachindex(m),
        kpi in eachindex(kp)
        
        if  col_int[kpi, mi] != 0
            wavespectrum[ii, jj, kk, kpi, mi] = was_copy[ii, jj, kk, kpi, mi] + dtau * col_int[kpi, mi]  
        end

    end

  
end