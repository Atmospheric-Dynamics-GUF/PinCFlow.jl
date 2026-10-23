"""
```julia 
StressTensor{A <: AbstractArray{<:AbstractFloat, 3}}
```

Container for the arrays to store the stress tensor elements.

```julia 
StressTensor(domain::Domain)::StressTensor
```

Construct a `StressTensor` instance with zero-initialized arrays.

# Fields

  - `pi11::A`: Stress tensor element ``\\Pi^{1 1}``.

  - `pi12::A`: Stress tensor element ``\\Pi^{1 2}``.
  
  - `pi13::A`: Stress tensor element ``\\Pi^{1 3}``.
  
  - `pi22::A`: Stress tensor element ``\\Pi^{2 2}``. 
  
  - `pi23::A`: Stress tensor element ``\\Pi^{2 3}``. 
  
  - `pi33::A`: Stress tensor element ``\\Pi^{3 3}``.

# Arguments

  - `domain`: Collection of domain-decomposition and MPI-communication parameters.
"""
struct StressTensor{A <: AbstractArray{<:AbstractFloat, 3}}
    pi11::A
    pi12::A
    pi13::A
    pi22::A
    pi23::A
    pi33::A
end

function StressTensor(domain::Domain)::StressTensor
    (; nxx, nyy, nzz) = domain

    return StressTensor((zeros(nxx, nyy, nzz) for i in 1:6)...)
end
