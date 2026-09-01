struct ChiParentTracker
    ntracked::Ref{Int}

    kpi_ref::Vector{Int}
    mi_ref::Vector{Int}

    kpi::Matrix{Int}
    mi::Matrix{Int}
    valid::BitMatrix

    kpi_radius::Vector{Int}
    mi_radius::Vector{Int}
end

function ChiParentTracker(wkb_mode::Union{NoWKB, SteadyState, SingleColumn, MultiColumn},
    triad_mode::NoTriad)::ChiParentTracker

    chi_tracker = ChiParentTracker(
        Ref(0),
        Int[],
        Int[],
        zeros(Int, 0, 0),
        zeros(Int, 0, 0),
        falses(0, 0),
        Int[],
        Int[],
    )
    return chi_tracker
end

function ChiParentTracker(namelists::Namelists,
    domain::Domain,
    wkb_mode::Union{NoWKB, SteadyState, SingleColumn, MultiColumn},
    triad_mode::Union{Triad2D, Triad3DIso})::ChiParentTracker

    (; wave_modes) = namelists.wkb
    (; nzz) = domain
    chi_tracker = ChiParentTracker(
    Ref(0),
    zeros(Int, wave_modes),
    zeros(Int, wave_modes),
    zeros(Int, wave_modes, nzz),
    zeros(Int, wave_modes, nzz),
    falses(wave_modes, nzz),
    zeros(Int, wave_modes),
    zeros(Int, wave_modes),
)
    return chi_tracker
end

