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
struct TriadNamelist{A <: Int, B <: Float64}
    kp_size::A
    m_size::A
    lkp::B
    lm::B
    triad_int::Bool
end

function TriadNamelist(;
    kp_size::Integer = 1,
    m_size::Integer = 1,
    lkp::Real = 1.0,
    lm::Real = 1.0, 
    triad_int = false,
)::TriadNamelist
    return TriadNamelist(
        Int(kp_size),
        Int(m_size),
        Float64(lkp),
        Float64(lm),
        Bool(triad_int)
    )
end
