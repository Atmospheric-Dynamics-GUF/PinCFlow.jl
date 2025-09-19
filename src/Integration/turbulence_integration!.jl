function turbulence_integration! end

function turbulence_integration!(
    state::State,
    p0::Predictands,
    dt::AbstractFloat,
)
    (; turbulencesetup) = state.namelists.turbulence

    turbulence_integration!(state, p0, dt, turbulencesetup)

    return
end

function turbulence_integration!(
    state::State,
    p0::Predictands,
    dt::AbstractFloat,
    turbulencesetup::NoTurbulence,
)
    return
end

function turbulence_integration!(
    state::State,
    p0::Predictands,
    dt::AbstractFloat,
    turbulencesetup::AbstractTurbulence,
)
    turbulence_integration!(state, p0, dt * 0.5, Dissipation())
    turbulence_integration!(state, p0, dt, Advection())
    turbulence_integration!(state, p0, dt, Diffusion())
    turbulence_integration!(state, p0, dt * 0.5, Dissipation())

    return
end

function turbulence_integration!(
    state::State,
    p0::Predictands,
    dt::AbstractFloat,
    process::Dissipation,
)
    (; alphaturb) = state.turbulence.turbulenceconstants
    (; tke) = state.turbulence.turbulencepredictands

    tke .= 4.0 ./ (alphaturb * dt .+ 2 ./ sqrt.(tke)) .^ 2
    return
end

function turbulence_integration!(
    state::State,
    p0::Predictands,
    dt::AbstractFloat,
    process::Advection,
)
    (; nstages, stepfrac) = state.time
    (; turbulencesetup) = state.namelists.turbulence

    for rkstage in 1:nstages
        reconstruct!(state)
        set_boundaries!(state, BoundaryReconstructions())

        compute_fluxes!(state, p0)

        set_boundaries!(state, BoundaryFluxes())

        save_backups!(state, :rho)

        update!(state, dt, rkstage, turbulencesetup)
        apply_unified_sponge!(
            state,
            stepfrac[rkstage] * dt,
            turbulencesetup,
        )
        set_boundaries!(state, BoundaryPredictands())
    end

    return

    # Advection and S+B-terms using RK3
    return
end

function turbulence_integration!(
    state::State,
    p0::Predictands,
    dt::AbstractFloat,
    process::Diffusion,
)

    (; nxx, nyy, nzz) = state.domain
    (; dz) = state.grid
    (; i0, i1, j0, j1, k0, k1) = state.domain
    (; k_diff) = state.turbulence.turbulenceconstants
    (; tke) = state.turbulence.turbulencepredictands
    (; dtke) = state.turbulence.turbulenceincrements
    
    # Define the intermediate variables used in the calculation

    ta = k1 - k0 + 1
    aa = zeros(Float64, ta)
    bb = zeros(Float64, ta)  
    cc = zeros(Float64, ta)  
    dd = zeros(Float64, ta)  
    qq = zeros(Float64, ta)
    ss = zeros(Float64, ta)
    K_ek = zeros(Float64, nxx, nyy, nzz)
    dr = dt / (2 * dz^2)

    # Defining the diffusion constant

    K_ek .= sqrt.(tke) .^ 3 .* k_diff

    # Defining the Thomas Algorithm for vertical diffusion 

    for j in j0:j1, i in i0:i1
                
        for k in k0:k1, 

            K_ekp = (K_ek[i, j, k] + K_ek[i, j, k+1])/2
            K_ekm = (K_ek[i, j, k-1] + K_ek[i, j, k])/2
                
            # Initialising the tridiagonal matrix

            aa[k-k0+1] = -dr * K_ekm
            bb[k-k0+1] = 1 + dr * K_ekm + dr * K_ekp
            cc[k-k0+1] = -dr * K_ekp

            # Initialising the RHS of the Algorithm

            dd[k-k0+1] = (1 - dr * K_ekm - dr * K_ekp) * tke[i, j, k] - (dr * K_ekp * tke[i, j, k+1]) - (dr * K_ekm * tke[i, j, k-1])

        end

        # Begin Thomas Algorithm

        dmx = dd[ta]
    
        # Forward Sweep

        qq[1] = -cc[1]/bb[1]
        dd[1] = dd[1]/bb[1]
        ss[1] = -aa[1]/bb[1]

        for mm in 2:ta
            qq[mm] = -cc[mm]/(bb[mm] + aa[mm]*qq[mm-1])
            dd[mm] = (dd[mm] - aa[mm]*dd[mm-1])/(bb[mm] + aa[mm]*qq[mm-1])
            ss[mm] = (-aa[mm] * ss[mm-1])/(bb[mm] + aa[mm]*qq[mm-1])
        end 

        # Backward Sweep

        qq[ta] = 0
        ss[ta] = 1 

        for mm in ta-1:-1:1
            ss[mm] = ss[mm] + qq[mm]*ss[mm+1]
            qq[mm] = dd[mm] + qq[mm]*qq[mm+1]
        end

        # Final Sweep

        dd[ta] = ( dmx - cc[ta]*qq[1] - aa[ta]*qq[ta-1] )/( cc[ta]*ss[1] + aa[ta]*ss[ta-1] + bb[ta] )

        for mm in 1:ta-1
            dd[mm] = dd[ta]*ss[mm] + qq[mm]
        end

        # End Thomas Algorithm

        # Setting the values obtained from the Thomas Algorithm to be the output values of this step
        for k in k0:k1, 
            tke[i, j, k] = dd[k-k0+1]
        end

    end

    return

end