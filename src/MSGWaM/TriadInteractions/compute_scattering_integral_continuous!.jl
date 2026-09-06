function compute_scattering_integral_continuous! end

function compute_scattering_integral_continuous!(
    state::State,
    ii::Integer,
    jj::Integer,
    kk::Integer,
    triad_mode::Triad2D,
)
    (; spec_tend) = state
    (; kp, m, kpl, kpc) = spec_tend.spec_grid
    (; aa, la, qq, lq, lia, liq, loglia, logliq) = spec_tend.kin_box
    (; wavespectrum, col_int, diag_dephasing_time) = spec_tend
    (; n2, rhobar) = state.atmosphere
    (; u) = state.variables.predictands
    (; compute_dephasing_time, nthreads_triad) = state.namelists.triad

    if compute_dephasing_time
        n_local = sqrt(n2[ii, jj, kk])
        dudz = compute_dphidz_center(u, state, ii, jj, kk, identity)
        dndz = compute_dphidz_center(n2, state, ii, jj, kk, sqrt)
    else
        n_local = 0.0
        dudz = 0.0
        dndz = 0.0
    end

    rhobar_local = rhobar[ii, jj, kk]

    was = @ivy view(wavespectrum, ii, jj, kk, :, :)

    update_interpolation_coef!(spec_tend, was, triad_mode)

    kpmin = kpc[1]
    kpmax = kpc[end]

    @sync for tid in 1:nthreads_triad
        inds = spec_tend.partition[tid]
        scr = spec_tend.scratch[tid]
        fpl, fpr, fq = scr.fpl, scr.fpr, scr.fq

        @spawn begin
            @ivy for idx in inds
                mi = (idx - 1) ÷ kpl + 1
                kpi = (idx - 1) % kpl + 1

                nk = was[kpi, mi]

                kr = kp[kpi]
                mr = m[mi]

                aar = aa[kpi]
                qqr = qq[kpi]

                fill!(view(fpl, 1:la[kpi]), 0.0)
                fill!(view(fpr, 1:la[kpi]), 0.0)
                fill!(view(fq, 1:lq[kpi]), 0.0)

                dephasing_time_k = Inf

                # ==========================================================
                # Sum interactions
                # ==========================================================

                sum_integral = 0.0

                if kr > 2.0 * kpmin
                    for i in 1:la[kpi]
                        pl = aar[i] - kr
                        pr = kr - aar[i]

                        st_value, tau_dep = compute_st_k(spec_tend, pl, 0.0, nk, kr, mr, n_local, dudz, dndz, 
                                                compute_dephasing_time, triad_mode, Sum())

                        fpl[i] = st_value
                        dephasing_time_k = min(dephasing_time_k, tau_dep)

                        if i == la[kpi]
                            fpr[i] = fpl[i]
                        else
                            st_value, tau_dep = compute_st_k(spec_tend, pr, 0.0, nk, kr, mr, n_local, dudz, dndz, 
                                                                compute_dephasing_time, triad_mode, Sum())

                            fpr[i] = st_value
                            dephasing_time_k = min(dephasing_time_k, tau_dep)
                        end
                    end

                    sum_integral =
                        trapazoidal_with_logbin(fpl, aar, la[kpi], lia[kpi], loglia[kpi]) +
                        trapazoidal_with_logbin(fpr, aar, la[kpi], lia[kpi], loglia[kpi])
                end

                # ==========================================================
                # Difference interactions
                # ==========================================================

                difference_integral = 0.0

                if kr < kpmax - kpmin
                    for j in 1:lq[kpi]
                        q = qqr[j]

                        st_value, tau_dep = compute_st_k(spec_tend, 0.0, q, nk, kr, mr, n_local, dudz, dndz, 
                                                        compute_dephasing_time, triad_mode, Difference())

                        fq[j] = st_value
                        dephasing_time_k = min(dephasing_time_k, tau_dep)
                    end

                    difference_integral = trapazoidal_with_logbin(fq, qqr, lq[kpi], liq[kpi], logliq[kpi])
                end

                col_int[ii, jj, kk, kpi, mi] =
                    2.0 * pi * (sum_integral - difference_integral) / rhobar_local

                diag_dephasing_time[kpi, mi] = dephasing_time_k
            end
        end
    end

    return nothing
end


function compute_scattering_integral!(
    state::State,
    ii::Integer,
    jj::Integer,
    kk::Integer,
    triad_mode::Triad3DIso
    )
    (; spec_tend) = state
    (; kp, m) = spec_tend.spec_grid
    (; aa, la, qq, lq, lia, liq, loglia, logliq) = spec_tend.kin_box
    (; wavespectrum, col_int) = spec_tend

    was = wavespectrum[ii, jj, kk, :, :]

    update_interpolation_coef!(spec_tend, was, triad_mode)

    for field in fieldnames(TriadTendencies)
        if field == :st_k || field == :col_int
            getfield(spec_tend, field) .= 0.0
        end
    end
    

     for mi in eachindex(m),
        kpi in eachindex(kp)

        nk = wavespectrum[ii, jj, kk, kpi, mi]

        kr = kp[kpi]  
        mr = m[mi]
        aar = aa[kpi]
        qqr = qq


        fpl = zeros(la[kpi])
        fql = zeros(lq)
        fpr = zeros(la[kpi])
        fqr = zeros(lq)

        chi = aar[1]/kr
        del = qqr[1]/kr

        for i in 1:(la[kpi]) #for p ∈ (-kr, kr)
            pl = aar[i] - kr  #for the left part of kinematic box, also as aar never equal to zero, so p = \pm 1 is not included here
            pr = kr - aar[i]  #for the right part of kinematic box

            for j in 1:lq #this excludes q = 0 as q_min \ne 0
                q = qqr[j]

                fql[j] = compute_st_k(spec_tend, pl, q, nk, kr, mr, triad_mode) / compute_delta_pq(kr, pl, q)
                fqr[j] = compute_st_k(spec_tend, pr, q, nk, kr, mr, triad_mode) / compute_delta_pq(kr, pr, q)
                
            end

            if i == la[kpi] # to avoid to count p = 0 twice
                fpl[i] = trapazoidal_with_logbin(fql, qqr, lq, liq, logliq) #integration of st_k w.r.t. q for fixed p=0
                fpr[i] = 0.0
            else
                fpl[i] = trapazoidal_with_logbin(fql, qqr, lq, liq, logliq) #integration of st_k w.r.t. q for fixed p in left plane
                fpr[i] = trapazoidal_with_logbin(fqr, qqr, lq, liq, logliq) #integration of st_k w.r.t. q for fixed p in right plane
            end


            # for the singularities at p = \pm kr, q = 0 as q = 0 was not included in the kinematic box

            fpl[i] +=  (compute_st_k(spec_tend, pl, 0.0, nk, kr, mr, triad_mode) + compute_st_k(spec_tend, pl, qqr[1], nk, kr, mr, triad_mode)) / sqrt(2*del/(kr^2-pl^2)) 

            fpr[i] +=  (compute_st_k(spec_tend, pr, 0.0, nk, kr, mr, triad_mode) + compute_st_k(spec_tend, pr, qqr[1], nk, kr, mr, triad_mode)) / sqrt(2*del/(kr^2-pr^2))
                

            
        end
        col_int[kpi, mi] = trapazoidal_with_logbin(fpl, aar, la[kpi], lia[kpi], loglia[kpi]) + 
                            trapazoidal_with_logbin(fpr, aar, la[kpi], lia[kpi], loglia[kpi])

        # Singularities p=±kr q≠0 
        for j in 1:lq
            q = qq[j]
            fql[j] = (compute_st_k(spec_tend, -kr, q, nk, kr, mr, triad_mode) + 
                compute_st_k(spec_tend, -kr * (1-chi), q, nk, kr, mr, triad_mode)+
                compute_st_k(spec_tend, kr * (1-chi), q, nk, kr, mr, triad_mode) +
                compute_st_k(spec_tend, kr, q, nk, kr, mr, triad_mode) ) *
                sqrt(2 * chi / (q * (2 * kr + q)))
        end

        col_int[kpi, mi] += trapazoidal_with_logbin(fql, qqr, lq, liq, logliq)

        # Singularities p=±kr q=0 here

        col_int[kpi, mi] += (compute_st_k(spec_tend, kr * (chi - 1), qq[1], nk, kr, mr, triad_mode) +
                            compute_st_k(spec_tend, kr * (1 - chi), qq[1], nk, kr, mr, triad_mode) ) *
                            2 * (pi - 2*asin(1-chi)) * asin(sqrt(del/2))

        
        col_int[kpi, mi] *= 4 * pi 

        
    end
   
end

