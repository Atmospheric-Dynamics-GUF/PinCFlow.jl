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
struct TriadNamelist{A <: Int, B <: Float64, C <: AbstractTriad, D <: AbstractTimeStepping}
    k_size::A
    l_size::A
    m_size::A
    lk::B
    ll::B
    lm::B
    triad_mode::C
    time_scheme::D
end

function TriadNamelist(;
    k_size::Integer = 1,
    l_size::Integer = 1,
    m_size::Integer = 1,
    lk::Real = 1.0,
    ll::Real = 1.0,
    lm::Real = 1.0, 
    triad_mode = NoTriad(),
    time_scheme = EulerMethod(),
)::TriadNamelist
    return TriadNamelist(
        Int(k_size),
        Int(l_size),
        Int(m_size),
        Float64(lk),
        Float64(ll),
        Float64(lm),
        triad_mode,
        time_scheme
    )
end
