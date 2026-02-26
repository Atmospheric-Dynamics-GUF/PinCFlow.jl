using Base.Threads: @spawn

function compute_scattering_integral! end

function compute_scattering_integral!(
    state::State,
    ii::Integer,
    jj::Integer,
    kk::Integer,
    triad_mode::Triad2D
    )
    (; spec_tend) = state.wkb
    (; kp, m, kpl, ml) = spec_tend.spec_grid
    (; aa, la, qq, lq, lia, liq, loglia, logliq) = spec_tend.kin_box
    (; wavespectrum, col_int) = spec_tend
    (; n2) = state.atmosphere

    nn = sqrt(n2[ii, jj, kk]) ##Bruint vaisala frequency at the level (ii, jj, kk)

    was = @ivy view(wavespectrum, ii, jj, kk, :, :)

    update_interpolation_coef!(spec_tend, was, triad_mode)

    ntasks = state.namelists.triad.nthreads_triad
    #@assert length(spec_tend.scratch) == ntasks
    #@assert length(spec_tend.partition) == ntasks

     

     @sync for tid in 1:ntasks
            inds = spec_tend.partition[tid]
            scr  = spec_tend.scratch[tid]
            fpl, fpr, fq = scr.fpl, scr.fpr, scr.fq

            @spawn begin
                @ivy for idx in inds
                    # map linear idx -> (kpi, mi)
                    mi  = (idx - 1) ÷ kpl + 1
                    kpi = (idx - 1) % kpl + 1



                    nk = was[kpi, mi]

                    kr = kp[kpi]  
                    mr = m[mi]
                    aar = aa[kpi]
                    qqr = qq[kpi]

                    # reset only the used parts
                    fill!(view(fpl, 1:la[kpi]), 0.0)
                    fill!(view(fpr, 1:la[kpi]), 0.0)
                    fill!(view(fq,  1:lq[kpi]), 0.0)

                    for i in 1:(la[kpi]) #for p ∈ (-kr, kr)
                        pl = aar[i] - kr  #for the left part of kinematic box, also as aar never equal to zero, so p = \pm 1 is not included here
                        pr = kr - aar[i]  #for the right part of kinematic box

                        if i == la[kpi] # to avoid to count p = 0 twice
                            fpl[i] = compute_st_k(spec_tend, pl, 0.0, nk, kr, mr, nn, triad_mode, Sum())
                            fpr[i] = 0.0
                        else
                            fpl[i] = compute_st_k(spec_tend, pl, 0.0, nk, kr, mr, nn, triad_mode, Sum())
                            fpr[i] = compute_st_k(spec_tend, pr, 0.0, nk, kr, mr, nn, triad_mode, Sum())
                        end    
                    end
                    #if kpi == kpl
                    #    col_int[ii, jj, kk, kpi, mi] = trapazoidal_with_logbin(fpl, aar, la[kpi], lia[kpi], loglia[kpi]) + 
                    #                                trapazoidal_with_logbin(fpr, aar, la[kpi], lia[kpi], loglia[kpi])
                    #else
                    for j in 1:(lq[kpi])
                        q = qqr[j]
                        fq[j] = compute_st_k(spec_tend, 0.0, q, nk, kr, mr, nn, triad_mode, Difference())
                    end
                    col_int[ii, jj, kk, kpi, mi] = trapazoidal_with_logbin(view(fpl,1:la[kpi]), aar, la[kpi], lia[kpi], loglia[kpi]) + 
                                                trapazoidal_with_logbin(view(fpr,1:la[kpi]), aar, la[kpi], lia[kpi], loglia[kpi]) -
                                                trapazoidal_with_logbin(view(fq, 1:lq[kpi]), qqr, lq[kpi], liq[kpi], logliq[kpi])
                    #end
                    # Singularities p=±kr, yet to define
                    col_int[ii, jj, kk, kpi, mi] *= (2 * pi) 
                end
            end
                
        end
   
end


function compute_scattering_integral!(
    state::State,
    ii::Integer,
    jj::Integer,
    kk::Integer,
    triad_mode::Triad3DIso
    )
    (; spec_tend) = state.wkb
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

