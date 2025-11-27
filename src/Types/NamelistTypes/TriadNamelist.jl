"""
```julia
TriadNamelist{A <: Int, B <: Float64, C <: MPI.Comm}
```

Namelist for parameters describing the triad interactions domain.

```julia
TriadNamelist(;
    kp_size::Integer = 1,
    m_size::Integer = 1,
    lkh::Real = 1,
    lm::Real = 1, 
)::TriadNamelist
```

Construct a `TriadNamelist` instance with the given keyword arguments as properties, converting them to meet the type constraints (Relevant only for the triad interactions only with WKB module).

# Fields/Keywords

  - `kp_size::A`: Number of grid cells in ``\\widehat{k_h}``-direction .

  - `kp_size::A`: Number of grid cells in ``\\widehat{m}``-direction (Relevant for the triad interactions only).

  - `lkh::B`: Domain extent in ``\\widehat{k_h}``-direction.

  - `lm::B`: Domain extent in ``\\widehat{m}``-direction.

"""
struct TriadNamelist{A <: Int, B <: Float64}
    kp_size::A
    m_size::A
    lkh::B
    lm::B
end

function TriadNamelist(;
    kp_size::Integer = 1,
    m_size::Integer = 1,
    lkh::Real = 1,
    lm::Real = 1, 
)::TriadNamelist
    return TriadNamelist(
        Int(kp_size),
        Int(m_size),
        Float64(lkh),
        Float64(lm)
    )
end
