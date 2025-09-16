"""
```julia
turbulence_computation!(
    state::State,
    p0::Predictands,
    dtstage::AbstractFloat,
    time::AbstractFloat,
    process::Dissipation,
)
```

Performs the dissipation step of the turbulence computation.

```julia
turbulence_computation!(
    state::State,
    p0::Predictands,
    dtstage::AbstractFloat,
    time::AbstractFloat,
    process::Diffusion,
)
```

Performs the diffusion step of the turbulence computation.

# Arguments

  - `state`: Model state.

  - `p0`: The turbulencepredictands that are used to compute the diffusion constant.

  - `dtstage`: Fractional time step.

  - `time`: Simulation time.

  - `process`: Denotes the process of the turbulence computation
"""
function turbulence_computation! end

function turbulence_computation!(
    state::State,
    p0::Predictands,
    dtstage::AbstractFloat,
    time::AbstractFloat,
    process::Dissipation,
)
    (; alphaturb) = state.turbulence.turbulenceconstants
    (; tke) = state.turbulence.turbulencepredictands
    (; dtke) = state.turbulence.turbulenceincrements
    (; nxx, nyy, nzz) = state.domain

    dis = zeros(Float64, nxx, nyy, nzz)
    
    if !(typeof(state.namelists.turbulence.turbulencesetup) <: NoTurbulence)
        tke .= 4.0 ./ (alphaturb * dtstage .+ 2 ./ sqrt.(tke)) .^ 2
    end

    return
end

function turbulence_computation!(
    state::State,
    p0::Predictands,
    dtstage::AbstractFloat,
    time::AbstractFloat,
    process::Diffusion,
)
    (; nxx, nyy, nzz) = state.domain
    (; dz) = state.grid
    (; i0, i1, j0, j1, k0, k1) = state.domain
    (; k_diff) = state.turbulence.turbulenceconstants
    (; tke) = state.turbulence.turbulencepredictands
    (; dtke) = state.turbulence.turbulenceincrements
    
    # p0_tke = p0.tke  tke old

    if !(typeof(state.namelists.turbulence.turbulencesetup) <: NoTurbulence)

        # Define the intermediate variables used in the calculation

        ta = k1 - k0 + 1
        aa = zeros(Float64, ta)
        bb = zeros(Float64, ta)  
        cc = zeros(Float64, ta)  
        dd = zeros(Float64, ta)  
        qq = zeros(Float64, ta)
        ss = zeros(Float64, ta)
        K_ek = zeros(Float64, nxx, nyy, nzz)
        dr = dtstage / (2 * dz^2)

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
            
    end

    return
end
