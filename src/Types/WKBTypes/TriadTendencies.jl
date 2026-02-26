"""
```julia
TriadTendencies{A <: AbstractVector{<:AbstractFloat}, B <: AbstractArray{<: AbstractFloat, 5}}
```

Triad interaction fields.

```julia
TriadTendencies{A <: AbstractVector{<:AbstractFloat}, B <: AbstractArray{<: AbstractFloat, 5}}
    kp::A
    m::A
    wavespectrum::B
    st_k::B
    col_int::B
    kin_box::KinematicBox
    interp_coef::InterpCoef
```

Construct a `TriadTendencies` instance, with arrays sized according to the given dimensions.

# Fields

  - `kp::A`: horizontal wave number grid.

  - `m::A`: vertical wave number grid.

  - `wavespectrum::B`: wave action spectrum in 5D phase space.
  
  - `st_k::B`: integrand of collision integral.

  - `col_int::B`: collision integral.

  - `kin_box::KinematicBox`: kinematic box container.


# Arguments

  - `domain`

  - `constants`
  
  - `wkb_mode`
  
"""

const RaySignature = Tuple{Int, Int, Int, Int, Float64, Float64, Float64}

struct TriadScratch{T}
    fpl::Vector{T}
    fpr::Vector{T}
    fq::Vector{T}
end

@inline function make_partition(N::Int, ntasks::Int)
    parts = Vector{UnitRange{Int}}(undef, ntasks)
    for tid in 1:ntasks
        a = (N * (tid - 1)) ÷ ntasks + 1
        b = (N * tid) ÷ ntasks
        parts[tid] = (a <= b) ? (a:b) : (1:0)   # empty chunk if needed
    end
    return parts
end


struct TriadTendencies{A <: AbstractArray{<: AbstractFloat, 5}, B <: AbstractArray{<: AbstractFloat, 2}}
    spec_grid::SpectralGrid
    wavespectrum::A
    ray_vol_signature::Array{Vector{RaySignature},5}
    col_int::A
    diag_time::B
    kin_box::KinematicBox
    interp_coef::InterpCoef
    res_manifold::ResManifold
    scratch::Vector{TriadScratch{Float64}}
    partition::Vector{UnitRange{Int}}
end


function TriadTendencies(namelists::Namelists,
   domain::Domain,
   constants::Constants,
   wkb_mode::Union{NoWKB, SteadyState, SingleColumn, MultiColumn},
   triad_mode::NoTriad)::TriadTendencies

    (; nthreads_triad) = namelists.triad

    spec_grid = SpectralGrid(wkb_mode, triad_mode)
    
    kin_box = KinematicBox(wkb_mode, triad_mode)

    interp_coef = InterpCoef(wkb_mode, triad_mode)
    
    res_manifold = ResManifold(wkb_mode, triad_mode)

    scratch = [TriadScratch(zeros(0), zeros(0), zeros(0))]

    partition = UnitRange{Int}[]
    
    return TriadTendencies(spec_grid, zeros(0, 0, 0, 0, 0), Array{Vector{RaySignature}}(undef, 0, 0, 0, 0, 0), zeros(0, 0, 0, 0, 0), zeros(0, 0),
      kin_box, interp_coef, res_manifold, scratch, partition)
end

function TriadTendencies(namelists::Namelists,
   domain::Domain,
   constants::Constants,
   wkb_mode::Union{SteadyState, SingleColumn, MultiColumn},
   triad_mode::Union{Triad2D, Triad3DIso})::TriadTendencies

    (; nxx, nyy, nzz) = domain
    (; rm_index, nthreads_triad) = namelists.triad

    # Compute the grid in kp and m direction

    spec_grid = SpectralGrid(namelists, constants, wkb_mode, triad_mode)
    
    (; kp, m, kpl, ml) = spec_grid

    # kinematic box and interpolation coefficients for 3D model
    #currently the model needs 
    if kpl > 1 && ml > 1
        kpmin = kp[1] 
        kpmax = kp[end]
          
        amin = (kpmin / kpl) .* ones(kpl)
        amax = Float64.(kp)
        ma = Int.(max.(8*ones(kpl), 1:kpl)) 
        if triad_mode == Triad3DIso() 
            mq = Int.(2 * kpl .* ones(kpl))
            qmin = ones(kpl) .* (kpmin / mq[1])
            qmax = ones(kpl) .* (2.0 * kpmax)
        else
            mq = reverse(Int.(max.(8*ones(kpl), 2 .* (1:kpl))))
            qmin = ones(kpl) .* (kpmin / mq[1])
            qmax = ones(kpl) .* (2.0 * kpmax)
        end

        kin_box = KinematicBox(amin, amax, ma, qmin, qmax, mq, wkb_mode, triad_mode)
      
    else
        error("Error in triad domain configurations, don't meet the specification for either 2D or 3D model")
    end

    res_manifold = ResManifold(kin_box.la[rm_index[1]], kin_box.lq[rm_index[1]], wkb_mode, triad_mode)

    interp_coef = InterpCoef(kp, m, wkb_mode, triad_mode)

    wavespectrum =  zeros(nxx, nyy, nzz, kpl, ml)
    col_int = zeros(nxx, nyy, nzz, kpl, ml)
    diag_time = zeros(kpl, ml)

    ray_vol_signature = Array{Vector{RaySignature},5}(undef, nxx, nyy, nzz, kpl, ml)

    for idx in eachindex(ray_vol_signature)
        ray_vol_signature[idx] = RaySignature[]
    end

    max_la = maximum(kin_box.la)
    max_lq = maximum(kin_box.lq)

    scratch = [TriadScratch(zeros(max_la), zeros(max_la), zeros(max_lq))
              for _ in 1:nthreads_triad]
    spec_l = kpl * ml
    partition = make_partition(spec_l, nthreads_triad) 

    return TriadTendencies(spec_grid, wavespectrum, ray_vol_signature, col_int, diag_time, kin_box, interp_coef, res_manifold, scratch, partition)

end

