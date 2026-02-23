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
struct TriadNamelist{A <: Int, B <: Float64, C <: AbstractTriad, D <: AbstractTimeStepping, E <: Tuple{Int, Int}}
    k_size::A
    l_size::A
    m_size::A
    k_max::B
    l_max::B
    m_max::B
    k_min::B
    l_min::B
    m_min::B
    col_int_tol::B
    triad_mode::C
    time_scheme::D
    rm_index::E
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
    col_int_tol::Real = 1.0E-5,
    triad_mode::AbstractTriad = NoTriad(),
    time_scheme::AbstractTimeStepping = EulerMethod(),
    rm_index::Tuple{Int, Int} = (1, 1),
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
        Float64(col_int_tol),
        triad_mode,
        time_scheme,
        rm_index
    )
end
