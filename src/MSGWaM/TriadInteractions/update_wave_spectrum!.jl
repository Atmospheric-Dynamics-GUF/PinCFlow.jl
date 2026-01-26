function update_wave_spectrum! end

function update_wave_spectrum!(
    state::State,
    ii::Integer,
    jj::Integer,
    kk::Integer,
    dtau::AbstractFloat,
    triad_mode::Union{Triad2D, Triad3DIso}
    )
    (; time_scheme) = state.namelists.triad
    update_wave_spectrum!(state, ii, jj, kk, dtau, triad_mode, time_scheme)
end


function update_wave_spectrum!(
    state::State,
    ii::Integer,
    jj::Integer,
    kk::Integer,
    dtau::AbstractFloat, 
    triad_mode::Union{Triad2D, Triad3DIso},
    time_scheme::EulerMethod,
)
    (; spec_tend) = state.wkb
    (; wavespectrum, col_int) = spec_tend
    (; kp, m) = spec_tend.spec_grid

    max_was = maximum(spec_tend.wavespectrum[ii, jj, kk, :, :])

    if max_was <= 1.0E-40
        #return if there is non significant wad in the physical grid cell
        return
    end

    compute_scattering_integral!(state, ii, jj, kk, triad_mode)
    tau_nl = get_nl_time_scale!(spec_tend, ii, jj, kk)

    if tau_nl > (dtau * 1.0E+5)
        #The nonlnear time scale is too large, inetraction is not required in this grid cell
        return
    end 


    #Euler method 
    @ivy for mi in eachindex(m),
        kpi in eachindex(kp)

        if  (dtau * abs(col_int[kpi, mi])) < (1.0E-5 * max_was) #cell wise exclusion for small col_int, the collision integral is too small just ignore it
            continue
        else
            wavespectrum[ii, jj, kk, kpi, mi] += dtau * col_int[kpi, mi]  
        end

    end
    
end

function update_wave_spectrum!(
    state::State,
    ii::Integer,
    jj::Integer,
    kk::Integer,
    dtau::AbstractFloat, 
    triad_mode::Union{Triad2D, Triad3DIso},
    time_scheme::Rk2Step,
)
    (; spec_tend) = state.wkb
    (; wavespectrum, col_int) = spec_tend
    (; kp, m) = spec_tend.spec_grid

    compute_scattering_integral!(state, ii, jj, kk, triad_mode)

    #RK2 method
    was_copy = spec_tend.wavespectrum
    @ivy for mi in eachindex(m),
        kpi in eachindex(kp)

        if  col_int[kpi, mi] != 0
            wavespectrum[ii, jj, kk, kpi, mi] += 0.5 * dtau * col_int[kpi, mi]
            #println("The collision integral is no zero for the grid cell", (ii, jj, kk, kpi, mi), 
            #"\n The collision integral is", col_int[kpi, mi])  
        end

    end

    compute_scattering_integral!(state, ii, jj, kk, triad_mode)

    @ivy for mi in eachindex(m),
        kpi in eachindex(kp)
        
        if  col_int[kpi, mi] != 0
            wavespectrum[ii, jj, kk, kpi, mi] = was_copy[ii, jj, kk, kpi, mi] + dtau * col_int[kpi, mi]  
        end

    end

  
end