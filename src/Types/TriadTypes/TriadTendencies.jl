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
        parts[tid] = (a <= b) ? (a:b) : (1:0)
    end

    return parts
end


struct TriadTendencies{
    A <: AbstractArray{<:AbstractFloat, 5},
    B <: AbstractArray{<:AbstractFloat, 2},
    C <: AbstractArray{Bool, 5},
    D <: AbstractArray{<:AbstractFloat, 3},
}
    spec_grid::SpectralGrid
    wavespectrum::A
    was_pred::A
    was_ray_signature::C
    col_int::A
    diag_time::B
    nl_time_scale::D
    dephasing_time::D
    chi_parent::Vector{Float64}
    chi_max::Ref{Float64}
    chi_tracker::ChiParentTracker
    kin_box::KinematicBox
    interp_coef::InterpCoef
    res_manifold::ResManifold
    scratch::Vector{TriadScratch{Float64}}
    partition::Vector{UnitRange{Int}}
    prev_dt::Ref{<:AbstractFloat}
    action_ref::Ref{<:AbstractFloat}
end

function TriadTendencies(namelists::Namelists,
    domain::Domain,
    constants::Constants,)::TriadTendencies
    (; wkb_mode) = namelists.wkb
    (; triad_mode) = namelists.triad

    return TriadTendencies(namelists, domain, constants, wkb_mode, triad_mode)
    
end

function TriadTendencies(
    namelists::Namelists,
    domain::Domain,
    constants::Constants,
    wkb_mode::Union{NoWKB, SteadyState, SingleColumn, MultiColumn},
    triad_mode::NoTriad,
)::TriadTendencies

    (; nthreads_triad) = namelists.triad

    spec_grid = SpectralGrid(wkb_mode, triad_mode)
    kin_box = KinematicBox(wkb_mode, triad_mode)
    interp_coef = InterpCoef(wkb_mode, triad_mode)
    res_manifold = ResManifold(wkb_mode, triad_mode)
    chi_tracker = ChiParentTracker(wkb_mode, triad_mode)

    scratch = [TriadScratch(zeros(0), zeros(0), zeros(0))]
    partition = UnitRange{Int}[]

    return TriadTendencies(
        spec_grid,
        zeros(0, 0, 0, 0, 0),
        zeros(0, 0, 0, 0, 0),
        falses(0, 0, 0, 0, 0),
        zeros(0, 0, 0, 0, 0),
        zeros(0, 0),
        zeros(0, 0, 0),
        zeros(0, 0, 0),
        zeros(0),
        Ref(0.0),
        chi_tracker,
        kin_box,
        interp_coef,
        res_manifold,
        scratch,
        partition,
        Ref(0.0),
        Ref(0.0),
    )
end


function TriadTendencies(
    namelists::Namelists,
    domain::Domain,
    constants::Constants,
    wkb_mode::Union{SteadyState, SingleColumn, MultiColumn},
    triad_mode::Union{Triad2D, Triad3DIso},
)::TriadTendencies

    (; nxx, nyy, nzz) = domain
    (; x_size) = namelists.domain
    (; rm_index, nthreads_triad) = namelists.triad
    (; wave_modes) = namelists.wkb

    # Compute the spectral grid.
    spec_grid = SpectralGrid(namelists, constants, wkb_mode, triad_mode)
    (; kp, m, kpc, kpl, ml) = spec_grid

    if kpl <= 1 || ml <= 1
        error(
            "Error in triad domain configurations, don't meet the specification " *
            "for either 2D or 3D model",
        )
    end

    # -------------------------------------------------------------------------
    # x_size = 1:
    #
    # Horizontal wavenumbers are discrete modes.
    # There is no continuous a/q integration in the k-direction.
    # -------------------------------------------------------------------------

    if x_size == 1

        if wkb_mode != SingleColumn()
            error("Triad interactions with x_size = 1 require SingleColumn WKB mode.")
        end

        kin_box = KinematicBox(wkb_mode, triad_mode)
        res_manifold = ResManifold(wkb_mode, triad_mode)

        scratch = [
            TriadScratch(zeros(0), zeros(0), zeros(0))
            for _ in 1:nthreads_triad
        ]

    # -------------------------------------------------------------------------
    # x_size > 1:
    #
    # Existing continuous/logarithmic k formulation.
    # This applies to both SingleColumn and MultiColumn.
    # -------------------------------------------------------------------------

    else

        amax = Float64.(kp)
        ma = Int.(max.(8 * ones(kpl), 1:kpl))

        if triad_mode == Triad3DIso()

            kpmin = kp[1]
            kpmax = kp[end]

            amin = (kpmin / (1.0 + eps())) .* ones(kpl)

            mq = Int.(2 * kpl .* ones(kpl))
            qmin = ones(kpl) .* (kpmin / mq[1])
            qmax = ones(kpl) .* (2.0 * kpmax)

        else

            # Use the actual outer cell edges as the represented
            # horizontal spectral domain.
            kpmin = kpc[1]
            kpmax = kpc[end]

            mq = reverse(Int.(max.(8 * ones(kpl), 2 .* (1:kpl))))

            amin = similar(kp)
            qmin = similar(kp)
            qmax = similar(kp)

            for kpi in eachindex(kp)
                kr = kp[kpi]

                # Sum interaction:
                #
                #     2 * kpmin <= a <= kr
                #
                # If kr <= 2 * kpmin, no resolved sum manifold exists.
                # A positive dummy interval is constructed here and is
                # skipped later in compute_scattering_integral!.
                amin[kpi] = kr > 2.0 * kpmin ? 2.0 * kpmin : 0.5 * kr

                # Difference interaction:
                #
                #     2 * kpmin <= q <= 2 * (kpmax - kr)
                #
                q_upper = 2.0 * (kpmax - kr)

                # If q_upper <= 2 * kpmin, no resolved difference
                # manifold exists. The dummy interval is skipped later.
                qmin[kpi] = q_upper > 2.0 * kpmin ? 2.0 * kpmin : 0.5 * q_upper
                qmax[kpi] = q_upper
            end
        end

        kin_box = KinematicBox(amin, amax, ma, qmin, qmax, mq, wkb_mode, triad_mode)

        res_manifold = ResManifold(
            kin_box.la[rm_index[1]],
            kin_box.lq[rm_index[1]],
            wkb_mode,
            triad_mode,
        )

        max_la = maximum(kin_box.la)
        max_lq = maximum(kin_box.lq)

        scratch = [
            TriadScratch(zeros(max_la), zeros(max_la), zeros(max_lq))
            for _ in 1:nthreads_triad
        ]
    end

    # Interpolation coefficients.
    interp_coef = InterpCoef(kp, m, wkb_mode, triad_mode)

    # Triad-interaction fields.
    wavespectrum = zeros(nxx, nyy, nzz, kpl, ml)
    was_pred = zeros(nxx, nyy, nzz, kpl, ml)
    was_ray_signature = falses(nxx, nyy, nzz, kpl, ml)
    col_int = zeros(nxx, nyy, nzz, kpl, ml)

    diag_time = zeros(kpl, ml)
    nl_time_scale = zeros(nxx, nyy, nzz)
    dephasing_time = zeros(nxx, nyy, nzz)

    chi_parent = fill(NaN, wave_modes)
    chi_max = Ref(NaN)
    chi_tracker = ChiParentTracker(namelists, domain, wkb_mode, triad_mode)

    # Thread partition of the complete spectral grid.
    spec_l = kpl * ml
    partition = make_partition(spec_l, nthreads_triad)

    prev_dt = Ref(0.0)
    action_ref = Ref(0.0)

    return TriadTendencies(
        spec_grid,
        wavespectrum,
        was_pred,
        was_ray_signature,
        col_int,
        diag_time,
        nl_time_scale,
        dephasing_time,
        chi_parent,
        chi_max,
        chi_tracker,
        kin_box,
        interp_coef,
        res_manifold,
        scratch,
        partition,
        prev_dt,
        action_ref,
    )
end