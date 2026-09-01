"""
```julia
TriadNamelist{A <: Int, B <: Float64, C <: MPI.Comm}
```

Namelist for parameters describing the triad interactions domain.

```julia
TriadNamelist(;
    kp_size::Integer = 1,
    m_size::Integer = 1,
    lkp::Real = 1,
    lm::Real = 1, 
)::TriadNamelist
```

Construct a `TriadNamelist` instance with the given keyword arguments as properties, converting them to meet the type constraints (Relevant only for the triad interactions only with WKB module).

# Fields/Keywords

  - `kp_size::A`: Number of grid cells in ``\\widehat{k_h}``-direction .

  - `m_size::A`: Number of grid cells in ``\\widehat{m}``-direction (Relevant for the triad interactions only).

  - `lkp::B`: Domain extent in ``\\widehat{k_perp}``-direction.

  - `lm::B`: Domain extent in ``\\widehat{m}``-direction.

"""
struct TriadNamelist{A <: Int, B <: Float64, C <: AbstractTriad, D <: AbstractTimeStepping, E <: Tuple{Int, Int}, F <: Bool, G <: AbstractMergeMode}
    k_size::A
    l_size::A
    m_size::A
    k_max::B
    l_max::B
    m_max::B
    k_min::B
    l_min::B
    m_min::B
    increment_rel_tol::B
    action_rel_tol::B
    triad_mode::C
    time_scheme::D
    rm_index::E
    nthreads_triad::A
    compute_dephasing_time::F
    smooth_wave_spectrum::F
    nl_cfl_number::B
    pl_cfl_number::B
    dt_growth_factor::B
    launch_rays_action_rel_tol::B
    discarded_action_fraction_tol::B
    projection_scheme::G
    chi_z_min::Float64
    chi_z_max::Float64
    chi_action_rel_tol::Float64
end

function TriadNamelist(;
    k_size::Integer = 1,
    l_size::Integer = 1,
    m_size::Integer = 1,
    k_max::Real = 1.0,
    l_max::Real = 1.0,
    m_max::Real = 1.0,
    k_min::Real = 1.0, #the dafault values of the minimum wave numbers is set in spectral grid type
    l_min::Real = 1.0,
    m_min::Real = 1.0, 
    increment_rel_tol::Real = 1.0E-4,
    action_rel_tol::Real = 1.0E-8,
    triad_mode::AbstractTriad = NoTriad(),
    time_scheme::AbstractTimeStepping = EulerMethod(),
    rm_index::Tuple{Int, Int} = (1, 1),
    nthreads_triad::Integer = 1,
    compute_dephasing_time::Bool = true,
    smooth_wave_spectrum::Bool = true,
    nl_cfl_number::Real = 5.0E-1,
    pl_cfl_number::Real = 5.0E-1,
    dt_growth_factor::Real = 1.25,
    launch_rays_action_rel_tol::Real = 1.0E-5,
    discarded_action_fraction_tol::Real = 1.0E-5,
    projection_scheme::AbstractMergeMode = ConstantWaveEnergy(),
    chi_z_min = 20.0E3,
    chi_z_max = 20.0E3,
    chi_action_rel_tol = 1.0E-2,
)::TriadNamelist
    return TriadNamelist(
        Int(k_size),
        Int(l_size),
        Int(m_size),
        Float64(k_max),
        Float64(l_max),
        Float64(m_max),
        Float64(k_min),
        Float64(l_min),
        Float64(m_min),
        Float64(increment_rel_tol),
        Float64(action_rel_tol),
        triad_mode,
        time_scheme,
        rm_index,
        Int(nthreads_triad),
        compute_dephasing_time,
        smooth_wave_spectrum,
        Float64(nl_cfl_number),
        Float64(pl_cfl_number),
        Float64(dt_growth_factor),
        Float64(launch_rays_action_rel_tol),
        Float64(discarded_action_fraction_tol),
        projection_scheme,
        chi_z_min,
        chi_z_max,
        chi_action_rel_tol,
    )
end
