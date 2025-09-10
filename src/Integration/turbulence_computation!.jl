"""
```julia
turbulence_computation!(
    state::State,
    p2::Predictands,
    dtstage::AbstractFloat,
    time::AbstractFloat,
    process::Dissipation,
)
```

Performs the dissipation step of the turbulence computation.

```julia
turbulence_computation!(
    state::State,
    p2::Predictands,
    dtstage::AbstractFloat,
    time::AbstractFloat,
    process::Diffusion,
)
```

Performs the diffusion step of the turbulence computation.

# Arguments

  - `state`: Model state.

  - `p0`: The predictands that are used to compute the transporting velocities in the computation of the fluxes.

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
        #dis .= (2 ./ sqrt.(tke)) .+ (alphaturb * dtstage)
        dis .= (2 ./ tke) .+ (alphaturb * dtstage)
        dis .= dis .* dis
        tke .= 4 ./ dis
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
    (; nzz) = state.domain
    (; i0, i1, j0, j1, k0, k1) = state.domain
    (; K_ek) = state.turbulence.turbulenceconstants
    (; tke) = state.turbulence.turbulencepredictands
    (; dtke) = state.turbulence.turbulenceincrements
    
    if !(typeof(state.namelists.turbulence.turbulencesetup) <: NoTurbulence)
        
        # Defining the Thomas Algorithm for vertical diffusion

            # Define the intermediate variables used in the calculation

                aa = zeros(Float64, nzz)
                bb = zeros(Float64, nzz)  
                cc = zeros(Float64, nzz)  
                dd = zeros(Float64, nzz)  
                qq = zeros(Float64, nzz)
                ss = zeros(Float64, nzz)  

                
            # Initialising the tridiagonal matrix
        
            for mm in 1:nzz
                aa[mm] = -dtstage * K_ek
                bb[mm] = 1 + 2 * dtstage * K_ek
                cc[mm] = -dtstage * K_ek
            end

            println(nzz)

            println(k0)

            println(k1)

            # Initialising the RHS of the Algorithm

            for j in j0:j1, i in i0:i1
                
                for k in k0:k1, 

                    dd[k] = tke[i, j, k] + (dtstage * K_ek * tke[i, j, k+1]) - (2 * dtstage * K_ek * tke[i, j, k]) + (dtstage * K_ek * tke[i, j, k-1])

                end

                    # Begin Thomas Algorithm

                    dmx = dd[nzz]
    
                    # Forward Sweep

                    qq[1] = -cc[1]/bb[1]
                    dd[1] = dd[1]/bb[1]
                    ss[1] = -aa[1]/bb[1]

                    for mm in 2:nzz
                        qq[mm] = -cc[mm]/(bb[mm] + aa[mm]*qq[mm-1])
                        dd[mm] = (dd[mm] - aa[mm]*dd[mm-1])/(bb[mm] + aa[mm]*qq[mm-1])
                        ss[mm] = (-aa[mm] * ss[mm-1])/(bb[mm] + aa[mm]*qq[mm-1])
                    end 

                    # Backward Sweep

                    qq[nzz] = 0
                    ss[nzz] = 1 

                    for mm in nzz-1:-1:1
                        ss[mm] = ss[mm] + qq[mm]*ss[mm+1]
                        qq[mm] = dd[mm] + qq[mm]*qq[mm+1]
                    end

                    # Final Sweep

                    dd[nzz] = ( dmx - cc[nzz]*qq[1] - aa[nzz]*qq[nzz-1] )/( cc[nzz]*ss[1] + aa[nzz]*ss[nzz-1] + bb[nzz] )

                    for mm in 1:nzz-1
                        dd[mm] = dd[nzz]*ss[mm] + qq[mm]
                    end

                    # End Thomas Algorithm

                    # Setting the values obtained from the Thomas Algorithm to be the output values of this step
                    tke[i, j, :] .= dd[:]

            



            end
            
            



    end

    return
end
